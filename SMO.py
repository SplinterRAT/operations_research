data = {"lambda": 0.95,  "mu": 1, "m": 0 }


def setRo(data):
    data["ro"] = data["lambda"] / data["mu"]


def setP0(data):
    ro = data["ro"]
    if ro == 1:
        result = (1 - ro) / (1 - ro ** (data["m"] + 2))
    else:
        result = 1 / (data["m"] + 2)
    data["p0"] = result


def getP(data, k):
    return k * data["p0"]

def setP(data, k):
    name = "p" + str(k)
    data[name] = getP(data, k)


def setPotkaz(data):
    data["potkaz"] = getP(data, 1)


def setQ(data):
    data["q"] = 1 - data["potkaz"]


def setPsystem(data):
    data["psystem"] = 1 - data["potkaz"]



def setA(data):
    data["A"] = data["lambda"] * data["q"]


def setNochered(data):
    m = data["m"]
    ro = data["ro"]
    p0 = data["p0"]
    if m == 1:
        chislitel =  ro * ro * (1 - (ro ** m) * (m + 1 - m * ro))
        znamenatel = (1 - ro ** (m - 2)) * (1 - ro)
        result = p0 * chislitel / znamenatel
    else:
        result = m * (m + 1) / (2 * m + 4)
    data["nochered"] = result


def setNobsl(data):
    data["nobsl"] = data["ro"] * data["q"]


def setNsystem(data):
    data["nsystem"] = data["nochered"] + data["nobsl"]


def setTochered(data):
    data["tochered"] = data["nochered"] / (data["lambda"] * data["psystem"])


def setTsystem(data):
    data["tsystem"] = data["nsystem"] / (data["lambda"] * data["psystem"])


def setAllParams(data):
    setRo(data)
    setP0(data)
    setPotkaz(data)
    setQ(data)
    setPsystem(data)
    setA(data)
    setNochered(data)
    setNobsl(data)
    setNsystem(data)
    setTochered(data)
    setTsystem(data)


def printParam(name, value):
    print("{} \t: {}".format(name, value))


def printAllParams(data):
    printParam("Коєфіціент використання обєкта        ", data["ro"])
    printParam("Вірогідність того, що лінія вільня     ", data["p0"])
    printParam("Вірогідність відмови в обслуговуванні        ", data["potkaz"])
    printParam("Вірогідність прийняття заявки в систему    ", data["q"])
    printParam("Відносна пропускна спроможність    ", data["psystem"])
    printParam("Абсолютна пропускна спроможність         ", data["A"])
    printParam("Среднє число заявок в черзі          ", data["nochered"])
    printParam("Среднє число заявок в обслуговуванні     ", data["nobsl"])
    printParam("Среднє число заявок в СМО               ", data["nsystem"])
    printParam("Средній час очікування заявки в черзі  ", data["tochered"])
    printParam("Средній час перебування заявки в системі", data["tsystem"])
setAllParams(data)
printAllParams(data)