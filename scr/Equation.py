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
    def __init__(self, tab):
        pass

    @abstractmethod
    def compute(self, tem):
        pass

    def getData(self, tem_room, nPoints, tem_convert):
        tem_min = self.tem_min
        tem_max = self.tem_max
        fid = self.fid

        if tem_convert != 0.0:
            tem_min = tem_min + tem_convert
            tem_max = tem_max + tem_convert

        X = getInterval(tem_room, nPoints, tem_min, tem_max)
        Y = [self.compute(i) for i in X]

        return X, Y, fid

    def getDataAt(self, tem):
        tem_min = self.tem_min
        tem_max = self.tem_max

        # all temperatures should be in degree Celsius
        tem = tem - 273.15

        if (tem < tem_min) or (tem > tem_max):
            return []

        val = self.compute(tem)
        return [val]


class nullEquation(Equation):
    def __init__(self, tab):
        pass

    def compute(self, tem):
        pass

    def getData(self, tem_room, nPoints, tem_convert):
        return [], [], []

    def getDataAt(self, tem):
        return []


class dnsEquation1(Equation):
    def __init__(self, tab):
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
    def __init__(self, tab):
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
    def __init__(self, tab):
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
    def __init__(self, tab):
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
    def __init__(self, tab):
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


class hcpEquation1(Equation):
    def __init__(self, tab):
        tem_cri = tab['tem_cri']
        a = tab['a_hcp_1']
        b = tab['b_hcp_1']
        c = tab['c_hcp_1']
        d = tab['d_hcp_1']
        tem_min = tab['tem_min_1']
        tem_max = tab['tem_max_1']

        columns = tab.columns
        if 'fid' in columns:
            fid = tab['fid'].values
        else:
            fid = ['']

        a = float(a)
        b = float(b)
        c = float(c)
        d = float(d)
        tem_min = float(tem_min)
        tem_max = float(tem_max)

        self.tem_cri = tem_cri
        self.fid = fid
        #
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.tem_min = tem_min
        self.tem_max = tem_max

    def getData(self, tem_room, nPoints, tem_convert):
        if self.a == '%':
            return [], [], []
        X, Y, fid = Equation.getData(self, tem_room, nPoints, tem_convert)
        return X, Y, fid

    def compute(self, tem):
        a = self.a
        b = self.b
        c = self.c
        d = self.d

        val = a
        val = val + b * math.pow(tem / 100, 1)
        val = val + c * math.pow(tem / 100, 2)
        val = val + d * math.pow(tem / 100, 3)
        val = val * 8.314472
        return val


class hcpEquation2(Equation):
    def __init__(self, tab):
        columns = tab.columns
        if 'fid' in columns:
            fid = tab['fid'].values
        else:
            fid = ['']

        a = tab['a_hcp_2'].values[0]
        b = tab['b_hcp_2'].values[0]
        c = tab['c_hcp_2'].values[0]
        d = tab['d_hcp_2'].values[0]
        tem_min = tab['tem_min_2'].values[0]
        tem_max = tab['tem_max_2'].values[0]

        if a != '%':
            a = float(a)
            b = float(b)
            c = float(c)
            d = float(d)
            tem_min = float(tem_min)
            tem_max = float(tem_max)

        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.tem_min = tem_min
        self.tem_max = tem_max
        self.fid = fid

    def getData(self, tem_room, nPoints, tem_convert):
        if self.a == '%':
            return [], [], []
        X, Y, fid = Equation.getData(self, tem_room, nPoints, tem_convert)
        return X, Y, fid

    def compute(self, tem):
        a = self.a
        b = self.b
        c = self.c
        d = self.d

        val = a
        val = val + b * math.pow(tem / 100, 1)
        val = val + c * math.pow(tem / 100, 2)
        val = val + d * math.pow(tem / 100, 3)
        val = val * 8.314472
        return val


class hcpEquation3(Equation):
    def __init__(self, tab):
        columns = tab.columns
        if 'fid' in columns:
            fid = tab['fid'].values
        else:
            fid = ['']

        a = tab['a_hcp_3'].values[0]
        b = tab['b_hcp_3'].values[0]
        c = tab['c_hcp_3'].values[0]
        d = tab['d_hcp_3'].values[0]
        e = tab['e_hcp_3'].values[0]
        f = tab['f_hcp_3'].values[0]
        tem_min = tab['tem_min_3'].values[0]
        tem_max = tab['tem_max_3'].values[0]
        tem_cri = tab['tem_cri'].values[0]

        if a != '%':
            a = float(a)
            b = float(b)
            c = float(c)
            d = float(d)
            e = float(e)
            f = float(f)
            tem_min = float(tem_min)
            tem_max = float(tem_max)
            tem_cri = float(tem_cri)

        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f
        self.tem_min = tem_min
        self.tem_max = tem_max
        self.tem_cri = tem_cri
        self.fid = fid

    def getData(self, tem_room, nPoints, tem_convert):
        if self.a == '%':
            return [], [], []
        X, Y, fid = Equation.getData(self, tem_room, nPoints, tem_convert)
        return X, Y, fid

    def compute(self, tem):
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        e = self.e
        f = self.f
        tem_cri = self.tem_cri

        tr = tem / tem_cri

        val = a * math.log(1 - tr)
        val = val + b / (1 - tr)
        val = val + c * math.pow(tr, 0)
        val = val + d * math.pow(tr, 1)
        val = val + e * math.pow(tr, 2)
        val = val + f * math.pow(tr, 3)
        val = val * 8.314472
        return val
