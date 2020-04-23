import numpy as np
from scipy.optimize import linprog
import nashpy as nash
from fractions import Fraction # so that numbers are not displayed in decimal.
from warnings import warn
from cvxopt.modeling import variable, op
import time
def w_scipy ():
    A = np.array([[-1, -1, -1], [-1,2, 0], [0, 0, -1], [-1, 0, 0], [0, -1, 0]])
    b = np.array([-1000, 0, -340, 0, 0])
    c = np.array([10,15,25])
    res = linprog(c,method = 'simplex', A_ub=A, b_ub=b,bounds=(0, None))
    print('Optimal value:', res.fun, '\nX:', res.x)
    return
def w_out_scipy():
    class Simplex(object):
        def __init__(self, num_vars, constraints, objective_function):
            """
            num_vars: Number of variables
            equations: A list of strings representing constraints
            each variable should be start with x followed by a underscore
            and a number
            eg of constraints
            ['1x_1 + 2x_2 >= 4', '2x_3 + 3x_1 <= 5', 'x_3 + 3x_2 = 6']
            Note that in x_num, num should not be more than num_vars.
            Also single spaces should be used in expressions.
            objective_function: should be a tuple with first element
            either 'min' or 'max', and second element be the equation
            eg 
            ('min', '2x_1 + 4x_3 + 5x_2')
            For solution finding algorithm uses two-phase simplex method
            """
            self.num_vars = num_vars
            self.constraints = constraints
            self.objective = objective_function[0]
            self.objective_function = objective_function[1]
            self.coeff_matrix, self.r_rows, self.num_s_vars, self.num_r_vars = self.construct_matrix_from_constraints()
            del self.constraints
            self.basic_vars = [0 for i in range(len(self.coeff_matrix))]
            self.phase1()
            r_index = self.num_r_vars + self.num_s_vars

            for i in self.basic_vars:
                if i > r_index:
                    raise ValueError("Infeasible solution")

            self.delete_r_vars()

            if 'min' in self.objective.lower():
                self.solution = self.objective_minimize()

            else:
                self.solution = self.objective_maximize()
            self.optimize_val = self.coeff_matrix[0][-1]

        def construct_matrix_from_constraints(self):
            num_s_vars = 0  # number of slack and surplus variables
            num_r_vars = 0  # number of additional variables to balance equality and less than equal to
            for expression in self.constraints:
                if '>=' in expression:
                    num_s_vars += 1

                elif '<=' in expression:
                    num_s_vars += 1
                    num_r_vars += 1

                elif '=' in expression:
                    num_r_vars += 1
            total_vars = self.num_vars + num_s_vars + num_r_vars

            coeff_matrix = [[Fraction("0/1") for i in range(total_vars+1)] for j in range(len(self.constraints)+1)]
            s_index = self.num_vars
            r_index = self.num_vars + num_s_vars
            r_rows = [] # stores the non -zero index of r
            for i in range(1, len(self.constraints)+1):
                constraint = self.constraints[i-1].split(' ')

                for j in range(len(constraint)):

                    if '_' in constraint[j]:
                        coeff, index = constraint[j].split('_')
                        if constraint[j-1] is '-':
                            coeff_matrix[i][int(index)-1] = Fraction("-" + coeff[:-1] + "/1")
                        else:
                            coeff_matrix[i][int(index)-1] = Fraction(coeff[:-1] + "/1")

                    elif constraint[j] == '<=':
                        coeff_matrix[i][s_index] = Fraction("1/1")  # add surplus variable
                        s_index += 1

                    elif constraint[j] == '>=':
                        coeff_matrix[i][s_index] = Fraction("-1/1")  # slack variable
                        coeff_matrix[i][r_index] = Fraction("1/1")   # r variable
                        s_index += 1
                        r_index += 1
                        r_rows.append(i)

                    elif constraint[j] == '=':
                        coeff_matrix[i][r_index] = Fraction("1/1")  # r variable
                        r_index += 1
                        r_rows.append(i)

                coeff_matrix[i][-1] = Fraction(constraint[-1] + "/1")

            return coeff_matrix, r_rows, num_s_vars, num_r_vars

        def phase1(self):
            # Objective function here is minimize r1+ r2 + r3 + ... + rn
            r_index = self.num_vars + self.num_s_vars
            for i in range(r_index, len(self.coeff_matrix[0])-1):
                self.coeff_matrix[0][i] = Fraction("-1/1")
            coeff_0 = 0
            for i in self.r_rows:
                self.coeff_matrix[0] = add_row(self.coeff_matrix[0], self.coeff_matrix[i])
                self.basic_vars[i] = r_index
                r_index += 1
            s_index = self.num_vars
            for i in range(1, len(self.basic_vars)):
                if self.basic_vars[i] == 0:
                    self.basic_vars[i] = s_index
                    s_index += 1

            # Run the simplex iterations
            key_column = max_index(self.coeff_matrix[0])
            condition = self.coeff_matrix[0][key_column] > 0

            while condition is True:

                key_row = self.find_key_row(key_column = key_column)
                self.basic_vars[key_row] = key_column
                pivot = self.coeff_matrix[key_row][key_column]
                self.normalize_to_pivot(key_row, pivot)
                self.make_key_column_zero(key_column, key_row)

                key_column = max_index(self.coeff_matrix[0])
                condition = self.coeff_matrix[0][key_column] > 0

        def find_key_row(self, key_column):
            min_val = float("inf")
            min_i = 0
            for i in range(1, len(self.coeff_matrix)):
                if self.coeff_matrix[i][key_column] > 0:
                    val = self.coeff_matrix[i][-1] / self.coeff_matrix[i][key_column]
                    if val <  min_val:
                        min_val = val
                        min_i = i
            if min_val == float("inf"):
                raise ValueError("Unbounded solution")
            if min_val == 0:
                warn("Dengeneracy")
            return min_i

        def normalize_to_pivot(self, key_row, pivot):
            for i in range(len(self.coeff_matrix[0])):
                self.coeff_matrix[key_row][i] /= pivot

        def make_key_column_zero(self, key_column, key_row):
            num_columns = len(self.coeff_matrix[0])
            for i in range(len(self.coeff_matrix)):
                if i != key_row:
                    factor = self.coeff_matrix[i][key_column]
                    for j in range(num_columns):
                        self.coeff_matrix[i][j] -= self.coeff_matrix[key_row][j] * factor

        def delete_r_vars(self):
            for i in range(len(self.coeff_matrix)):
                non_r_length = self.num_vars + self.num_s_vars + 1
                length = len(self.coeff_matrix[i])
                while length != non_r_length:
                    del self.coeff_matrix[i][non_r_length-1]
                    length -= 1

        def update_objective_function(self):
            objective_function_coeffs = self.objective_function.split()
            for i in range(len(objective_function_coeffs)):
                if '_' in objective_function_coeffs[i]:
                    coeff, index = objective_function_coeffs[i].split('_')
                    if objective_function_coeffs[i-1] is '-':
                        self.coeff_matrix[0][int(index)-1] = Fraction(coeff[:-1] + "/1")
                    else:
                        self.coeff_matrix[0][int(index)-1] = Fraction("-" + coeff[:-1] + "/1")

        def check_alternate_solution(self):
            for i in range(len(self.coeff_matrix[0])):
                if self.coeff_matrix[0][i] and i not in self.basic_vars[1:]:
                    warn("Alternate Solution exists")
                    break

        def objective_minimize(self):
            self.update_objective_function()

            for row, column in enumerate(self.basic_vars[1:]):
                if self.coeff_matrix[0][column] != 0:
                    self.coeff_matrix[0] = add_row(self.coeff_matrix[0], multiply_const_row(-self.coeff_matrix[0][column], self.coeff_matrix[row+1]))

            key_column = max_index(self.coeff_matrix[0])
            condition = self.coeff_matrix[0][key_column] > 0

            while condition is True:

                key_row = self.find_key_row(key_column = key_column)
                self.basic_vars[key_row] = key_column
                pivot = self.coeff_matrix[key_row][key_column]
                self.normalize_to_pivot(key_row, pivot)
                self.make_key_column_zero(key_column, key_row)

                key_column = max_index(self.coeff_matrix[0])
                condition = self.coeff_matrix[0][key_column] > 0

            solution = {}
            for i, var in enumerate(self.basic_vars[1:]):
                if var < self.num_vars:
                    solution['x_'+str(var+1)] = self.coeff_matrix[i+1][-1]

            for i in range(0, self.num_vars):
                if i not in self.basic_vars[1:]:
                    solution['x_'+str(i+1)] = Fraction("0/1")
            self.check_alternate_solution()
            return solution

        def objective_maximize(self):
            self.update_objective_function()

            for row, column in enumerate(self.basic_vars[1:]):
                if self.coeff_matrix[0][column] != 0:
                    self.coeff_matrix[0] = add_row(self.coeff_matrix[0], multiply_const_row(-self.coeff_matrix[0][column], self.coeff_matrix[row+1]))

            key_column = min_index(self.coeff_matrix[0])
            condition = self.coeff_matrix[0][key_column] < 0

            while condition is True:

                key_row = self.find_key_row(key_column = key_column)
                self.basic_vars[key_row] = key_column
                pivot = self.coeff_matrix[key_row][key_column]
                self.normalize_to_pivot(key_row, pivot)
                self.make_key_column_zero(key_column, key_row)

                key_column = min_index(self.coeff_matrix[0])
                condition = self.coeff_matrix[0][key_column] < 0

            solution = {}
            for i, var in enumerate(self.basic_vars[1:]):
                if var < self.num_vars:
                    solution['x_'+str(var+1)] = self.coeff_matrix[i+1][-1]

            for i in range(0, self.num_vars):
                if i not in self.basic_vars[1:]:
                    solution['x_'+str(i+1)] = Fraction("0/1")

            self.check_alternate_solution()

            return solution

    def add_row(row1, row2):
        row_sum = [0 for i in range(len(row1))]
        for i in range(len(row1)):
            row_sum[i] = row1[i] + row2[i]
        return row_sum

    def max_index(row):
        max_i = 0
        for i in range(0, len(row)-1):
            if row[i] > row[max_i]:
                max_i = i

        return max_i

    def multiply_const_row(const, row):
        mul_row = []
        for i in row:
            mul_row.append(const*i)
        return mul_row

    def min_index(row):
        min_i = 0
        for i in range(0, len(row)):
            if row[min_i] > row[i]:
                min_i = i

        return min_i
    objective = ('maximize', '3x_1 + 4x_2')
    constraints = ['4x_1 + 1x_2 <= 8', '1x_1 - 1x_2 >= -3']
    Lp_system = Simplex(num_vars=2, constraints=constraints, objective_function=objective)
    print(Lp_system.solution)
    print(Lp_system.optimize_val)
def transport_problem():
    start = time.time()
    x = variable(9, 'x')
    c= [7,3,6,4,8,2,1,5,9]
    z=(c[0]*x[0] + c[1]*x[1] +c[2]* x[2] +c[3]*x[3] + c[4]*x[4] +c[5]* x[5]+c[6]*x[6] +c[7]*x[7] +c[8]* x[8])
    mass1 = (x[0] + x[1] +x[2] <= 74)
    mass2 = (x[3] + x[4] +x[5] <= 40)
    mass3 = (x[6] + x[7] + x[8] <= 36)
    mass4 = (x[0] + x[3] + x[6] == 20)
    mass5 = (x[1] +x[4] + x[7] == 45)
    mass6 = (x[2] + x[5] + x[8] == 30)
    mass7 = (x[1] == 30)
    x_non_negative = (x >= 0)    
    problem =op(z,[mass1,mass2,mass3,mass4 ,mass5,mass6, mass7,x_non_negative])
    problem.solve(solver='glpk')  
    print("Result:")
    for i in x.value:
            print(i)
    print("Costs:")
    print(problem.objective.value()[0])
    stop = time.time()
    print ("Time :")
    print(stop - start)
    return
def game_theory():
    A = np.array([[0, -1, 1], [1, 0, -1], [-1, 1, 0]])
    B = - A
    rps = nash.Game(A, B)
    print(rps)
    #utility of a pair strategies
    sigma_c = [1 / 2, 1 / 2, 0]
    sigma_r = [0, 1 / 2, 1 / 2]
    print(rps[sigma_r, sigma_c])
    #Nash equilibria
    eqs = rps.support_enumeration()
    print(list(eqs))
    return
transport_problem()
