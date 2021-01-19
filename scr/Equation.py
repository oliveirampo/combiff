from abc import ABC, abstractmethod
import numpy as np
import math


def getInterval(tem_room, n, tem_min, tem_max):
    rng = np.linspace(tem_min, tem_max, n)
    if tem_room > tem_min and tem_room < tem_max:
        idx = (np.abs(rng - tem_room)).argmin()
        rng[idx] = tem_room
    return rng


class Equation(ABC):
    @abstractmethod
    def compute(self, tem):
        pass


    @abstractmethod
    def addCoefficients(self, tab):
        pass


    def getData(self, tem_room, nPoints, tab, tem_convert):
        self.addCoefficients(tab)
        tem_min = self.tem_min
        tem_max = self.tem_max
        fid = self.fid

        if tem_convert != 0.0:
            tem_min = tem_min + tem_convert
            tem_max = tem_max + tem_convert

        X = getInterval(tem_room, nPoints, tem_min, tem_max)
        Y = [self.compute(i) for i in X]

        return X, Y, fid


    def getDataAt(self, tem, tab):
        self.addCoefficients(tab)
        tem_min = self.tem_min
        tem_max = self.tem_max

        # all temperaturess should be in degree Celsius
        tem = tem - 273.15

        if (tem < tem_min) or (tem > tem_max):
            return []

        val = self.compute(tem)
        return [val]


class nullEquation(Equation):
    def __init__(self):
        pass


    def compute(self, tem):
        pass


    def addCoefficients(self, tab):
        pass


    def getData(self, tem_room, nPoints, tab, tem_convert):
        return [], [], []


    def getDataAt(self, tem, tab):
        return []



class dnsEquation1(Equation):
    def __init__(self):
        pass

    def addCoefficients(self, tab):
        a = tab['a_dns']
        b = tab['b_dns']
        c = tab['c_dns']
        n = tab['n_dns']
        tem_min = tab['tem_min']
        tem_max = tab['tem_max']
        fid = tab['fid'].values

        a = float(a)
        b = float(b)
        c = float(c)
        n = float(n)
        tem_min = float(tem_min)
        tem_max = float(tem_max)

        self.a = a
        self.b = b
        self.c = c
        self.n = n
        self.tem_min = tem_min
        self.tem_max = tem_max
        self.fid = fid


    def compute(self, tem):
        a = self.a
        b = self.b
        c = self.c
        n = self.n
        exp = -math.pow(1.0 - (tem / c), n)
        return 1000.0 * a * math.pow(b, exp)


class hvpEquation1(Equation):
    def __init__(self):
        pass


    def addCoefficients(self, tab):
        a = tab['a_hvp']
        n = tab['n_hvp']
        tem_cri = tab['tem_cri']
        tem_min = tab['tem_min']
        tem_max = tab['tem_max']

        columns = tab.columns
        if 'fid' in columns:
            fid = tab['fid'].values
        else:
            fid = tab['met'].values

        a = float(a)
        n = float(n)
        tem_cri = float(tem_cri)
        tem_min = float(tem_min)
        tem_max = float(tem_max)

        self.a = a
        self.n = n
        self.tem_cri = tem_cri
        self.tem_min = tem_min
        self.tem_max = tem_max
        self.fid = fid


    def compute(self, tem):
        a = self.a
        n = self.n
        tem_cri = self.tem_cri
        return a * math.pow(1.0 - (tem / tem_cri), n)


class pvpEquation1(Equation):
    def __init__(self):
        pass


    def addCoefficients(self, tab):
        a = tab['a_pvp']
        b = tab['b_pvp']
        c = tab['c_pvp']
        tem_min = tab['tem_min']
        tem_max = tab['tem_max']

        columns = tab.columns
        if 'fid' in columns:
            fid = tab['fid'].values
        else:
            fid = tab['met'].values

        a = float(a)
        b = float(b)
        c = float(c)
        tem_min = float(tem_min)
        tem_max = float(tem_max)

        self.a = a
        self.b = b
        self.c = c
        self.tem_min = tem_min
        self.tem_max = tem_max
        self.fid = fid


    def compute(self, tem):
        # Torr
        a = self.a
        b = self.b
        c = self.c
        logP = a - ( b / (tem + c) )

        # Bar
        P = math.pow(10, logP) * 1.01325/760.0

        # kPA
        P = P * 100
        return P


class gamEquation1(Equation):
    def __init__(self):
        pass

    def addCoefficients(self, tab):
        a = tab['a_gam']
        b = tab['b_gam']
        n = tab['n_gam']
        tem_min = tab['tem_min']
        tem_max = tab['tem_max']

        columns = tab.columns
        if 'fid' in columns:
            fid = tab['fid'].values
        else:
            fid = tab['met'].values

        a = float(a)
        b = float(b)
        n = float(n)
        tem_min = float(tem_min)
        tem_max = float(tem_max)

        self.a = a
        self.b = b
        self.n = n
        self.tem_min = tem_min
        self.tem_max = tem_max
        self.fid = fid

    def compute(self, tem):
        a = self.a
        b = self.b
        n = self.n
        return a * (1 - tem / b) ** n


class alpEquation1(Equation):
    def __init__(self):
        pass

    def addCoefficients(self, tab):
        a = tab['a_alp']
        b = tab['b_alp']
        n = tab['n_alp']
        tem_min = tab['tem_min']
        tem_max = tab['tem_max']

        columns = tab.columns
        if 'fid' in columns:
            fid = tab['fid'].values
        else:
            fid = tab['met'].values

        a = float(a)
        b = float(b)
        n = float(n)
        tem_min = float(tem_min)
        tem_max = float(tem_max)

        self.a = a
        self.b = b
        self.n = n
        self.tem_min = tem_min
        self.tem_max = tem_max
        self.fid = fid

    def compute(self, tem):
        a = self.a
        b = self.b
        n = self.n
        exp = math.pow(1 - (tem / b), n)
        val = a * exp
        return val
