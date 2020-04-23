import math
import geojsonio
import matplotlib.pyplot as plt
import matplotlib
import shapely
from shapely.geometry import MultiPoint, Polygon, LineString , Point
import geopandas as gpd
import random
import itertools
from matplotlib.patches import Rectangle
import numpy as np
def start ():
    print ("Choose what square you want to calculate \n 1 - Function \n 2 - Country")
    chck = int(input())
    if chck == 1 :
        square_figure()
    elif chck == 2 :
        country_info()
    else :
        print ("Error")
#init for country
def init_country():
    print ("Enter the country name")
    country_name = input()
    print ("Enter number of dots")
    n_dots = int(input())
    return country_name, n_dots
#init for figure
def init_figure():
     #init
    print ("Enter lower bound")
    l_b = int(input())
    print ("Enter higher bound")
    h_b = int(input())
    print ("Accuracy of integral")
    N = int(input())
    print ("Enter number of dots")
    n_dots = int(input())
    return l_b, h_b, N , n_dots
#for points
def country_info ():
    country_name , n_dots = init_country()
    #reading file & getting country data
    countries = gpd.read_file('/home/splinter/Documents/Operations_Research/countries.geojson')
    country_data = countries[countries.ADMIN == country_name]
    #for reusing dots-function
    figure = country_data
    #parametrs of rectangle
    minmax = country_data.total_bounds
    minx , miny , maxx , maxy = minmax[0] , minmax[1] , minmax[2] , minmax[3]
    param_width  = math.sqrt(((maxx-minx)**2 + (miny-miny)**2))
    param_height = math.sqrt(((minx-minx)**2 + (maxy-miny)**2))
    #creating a rectangle
    rect = plt.Rectangle((minx , miny), param_width, param_height, edgecolor='r', facecolor="none")
    #creating plot
    fig = plt.figure()
    ax = fig.add_subplot()
    #creating plot for visualisation
    reslt = plt.figure()
    bx = reslt.add_subplot()
    #creating dots
    total_p = 0
    fig_p = 0
    rect_p = 0
    for i in range(n_dots) :
        x_dot = random.uniform(minx,maxx)
        y_dot = random.uniform(miny,maxy)
        c_in = figure.contains(Point(x_dot , y_dot)).iloc[0]
        #c_in_fil = list(filter(bool, c_in))
        r_in = True if minx <= x_dot <= maxx and miny <= y_dot <=maxy else False
        #print ("r_in " + str(r_in) + " c_in " + str(c_in))
        #test if it works
        if r_in == True and c_in == False:
            rect_p += 1
        elif r_in == True and c_in == True:
            fig_p += 1
        #drawing dots
        ax.plot(x_dot,y_dot, 'o')
    #Monte_Carlo
    print (rect_p)
    print (fig_p)
    total_p = rect_p + fig_p
    square = (fig_p / total_p) * (param_width * param_height)
    #square_real = square / 0.000187638
    plt.title(country_name, fontsize=19)
    plt.suptitle('S = ' + str(square), fontsize=12)
    #drawing shape
    ax.add_patch(rect)
    country_data.boundary.plot(ax=ax)
    country_data.boundary.plot(ax=bx)
    plt.show()
    #print(country_data)
    #print(minmax)
    print(param_width)
#1st part
def square_figure():
    l_b , h_b , N , n_dots = init_figure()
    #creating sinusoid
    x = np.arange(l_b,h_b, 0.1)
    y = np.sin(x) + 1
    #creating plot
    sin_pl = plt.figure()
    cx = sin_pl.add_subplot()
    #another plot for result
    reslt = plt.figure()
    dx = reslt.add_subplot()
    #rectangle
    minx = min(x)
    maxx = max(x)
    miny = min(y)
    maxy = max(y)
    param_width  = math.sqrt(((maxx-minx)**2 + (miny-miny)**2))
    param_height = math.sqrt(((minx-minx)**2 + (maxy-miny)**2))
    rect = plt.Rectangle((minx , miny), param_width, param_height, edgecolor='r', facecolor="none")
    #creating fiigure from plot
    verts = [(minx, 0), *zip(x, y), (maxx, 0)]
    poly = plt.Polygon(verts, facecolor='0.9', edgecolor='0.5')
    poly_check = Polygon(verts)
    dx.add_patch(poly)
    # !NEEED 1 more patch 4 cx
    #Monte Carlo
    square = 0
    total_p = 0
    fig_p = 0
    rect_p = 0
    for i in range(n_dots) :
        x_dot = random.uniform(minx,maxx)
        y_dot = random.uniform(miny,maxy)
        c_in = poly_check.contains(Point(x_dot , y_dot))
        #c_in_fil = list(filter(bool, c_in))
        r_in = True if minx <= x_dot <= maxx and miny <= y_dot <=maxy else False
        #print ("r_in " + str(r_in) + " c_in " + str(c_in))
        #test if it works
        if r_in == True and c_in == False:
            rect_p += 1
        elif r_in == True and c_in == True:
            fig_p += 1
        #drawing dots
        dx.plot(x_dot,y_dot, 'o')
    total_p = rect_p + fig_p
    square = (fig_p / total_p) * (param_width * param_height)
    print (total_p)
    print (fig_p)
    #Integral
    def integrate(f, a, b, N):
        x = np.linspace(a+(b-a)/(2*N), b-(b-a)/(2*N), N)
        fx = f (x) + 1
        area = np.sum(fx)*(b-a)/N
        return area
    square_i = integrate(np.sin, l_b, h_b, N)
    #Monte Carlo
    #drawing
    dx.add_patch(rect)
    dx.plot(x,y)
    cx.plot(x,y)
    reslt.suptitle("f(x) = sin(x)", fontsize=10)
    dx.set_title('S (Monte Carlo) = ' + str(square) + '\n' + 'S (integtal) = ' + str(square_i), fontsize=10)
    sin_pl.suptitle("f(x) = sin(x)", fontsize=10)
    cx.set_title('S (Monte Carlo) = ' + str(square) + '\n' + 'S (integtal) = ' + str(square_i), fontsize=10)
    plt.show()
start()