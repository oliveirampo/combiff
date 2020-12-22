from abc import ABC, abstractmethod
import numpy as np
import math
import sys


class Equation(ABC):
    @abstractmethod
    def compute(self, tem):
        pass

    @abstractmethod
    def getData(self, tem_room, nPoints, tab, tem_convert):
        pass


class nullEquation(Equation):
    def __init__(self):
        pass


    def compute(self, tem):
        pass


    def getData(self, tem_room, nPoints, tab, tem_convert):
        return [], [], []


class dnsEquation1(Equation):
    def __init__(self):
        pass

    def addCoefficients(self, a, b, c, n):
        self.a = a
        self.b = b
        self.c = c
        self.n = n


    def compute(self, tem):
        a = self.a
        b = self.b
        c = self.c
        n = self.n
        exp = -math.pow(1.0 - (tem / c), n)
        return 1000.0 * a * math.pow(b, exp)


    def getData(self, tem_room, nPoints, tab, tem_convert):
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

        self.addCoefficients(a, b, c, n)

        if tem_convert != 0.0:
            tem_min = tem_min + tem_convert
            tem_max = tem_max + tem_convert

        X = getInterval(tem_room, nPoints, tem_min, tem_max)
        Y = [self.compute(i) for i in X]

        return X, Y, fid


class hvpEquation1(Equation):
    def __init__(self):
        pass


    def addCoefficients(self, a, n, tem_cri):
        self.a = a
        self.n = n
        self.tem_cri = tem_cri


    def compute(self, tem):
        a = self.a
        n = self.n
        tem_cri = self.tem_cri
        return a * math.pow(1.0 - (tem / tem_cri), n)


    def getData(self, tem_room, nPoints, tab, tem_convert):
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

        self.addCoefficients(a, n, tem_cri)

        if tem_convert != 0.0:
            tem_min = tem_min + tem_convert
            tem_max = tem_max + tem_convert

        X = getInterval(tem_room, nPoints, tem_min, tem_max)
        Y = [self.compute(i) for i in X]

        return X, Y, fid


def getInterval(tem_room, n, tem_min, tem_max):
    rng = np.linspace(tem_min, tem_max, n)
    if tem_room > tem_min and tem_room < tem_max:
        idx = (np.abs(rng - tem_room)).argmin()
        rng[idx] = tem_room
    return rng

