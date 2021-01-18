from abc import ABC, abstractmethod
import numpy as np
import sys


import myExceptions


# if pre == '%' -> pre = 1.0 kPa
def checkPressure(col, data, defaultPressure):
    for i in range(data[:, col].shape[0]):
        try:
            data[i, col] = float(data[i, col])

        except ValueError:
            val = data[i, col]
            if val == '%':
                data[i, col] = defaultPressure
            else:
                raise myExceptions.WrongProperty('pressure', val, 'dbsRelation::checkPressure()')

    # data[:,col] = data[:,col].astype(np.float)
    return data


# if tem == '%' -> remove data
def checkTemperature(col, data):
    rm = []

    for i in range(data[:, col].shape[0]):
        try:
            data[i, col] = float(data[i, col])

        except ValueError:
            val = data[i, col]
            if val == 'sat._liquid':
                rm.append(i)
            else:
                raise myExceptions.WrongProperty('temperature', val, 'dbsRelation:checkTemperaure()')

    for idx in reversed(rm):
        data = np.delete(data, idx, 0)

    return data


# if val == '%' -> remove data
def checkProperty(col, data, propVar):
    for i in range(data[:, col].shape[0]):
        try:
            data[i, col] = float(data[i, col])

        except ValueError:
            val = data[i, col]
            if '+/-' in val:
                val = val.split('+/-')
                val = float(val[0])
                data[i, col] = val
            else:
                print(data[i, col])
                raise myExceptions.WrongProperty(propVar, val, 'dbsRelation::checkProperty()')


class DbsRelation(ABC):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):
        self.larsCode = larsCode
        self.prop = prop
        self.rel = rel
        self.var = var
        self.col = col
        #
        self.pre = pre
        self.tem = tem
        #
        self.fid = fid
        self.met = met
        #
        self.prop_convert = prop_convert
        self.tem_convert = tem_convert
        self.pre_convert = pre_convert
        #
        self.marker = marker
        self.color = color
        #
        self.preCol = 0
        self.temCol = 1
        self.propCol = 2
        self.fidCol = 3
        self.metCol = 4

    @abstractmethod
    def getLabel(self):
        pass

    def getValues(self, data):
        pre = data[:, self.preCol]
        tem = data[:, self.temCol]
        val = data[:, self.propCol]
        fid = data[:, self.fidCol]
        met = data[:, self.metCol]

        pre = pre.astype(np.float)
        tem = tem.astype(np.float)
        val = val.astype(np.float)

        return pre, tem, val, fid, met, self.marker, self.color

    def getData(self, tab, defaultPressure):
        data = self.getDataHelper(tab, defaultPressure)
        data = self.checkData(data, defaultPressure)

        col = self.temCol
        self.convertTemperature(col, data)

        col = self.propCol
        self.convertProperty(col, data)

        return data

    def getDataHelper(self, tab, defaultPressure):
        columns = self.col.strip().split()

        if columns[1] != 'rec':
            print('Problem with col sequence for {}_{}'.format(self.rel, self.larsCode))
            sys.exit(666)

        preVar = self.pre
        temVar = self.tem
        propVar = self.var
        fidVar = self.fid
        metVar = self.met

        prop = tab[propVar].to_numpy()

        pre = prop.shape[0] * [defaultPressure]
        if preVar != '':
            pre = tab[preVar].to_numpy()

        tem = prop.shape[0] * [298.15]
        if temVar != '':
            tem = tab[temVar].to_numpy()

        fid = prop.shape[0] * ['']
        if fidVar != '':
            fid = tab[fidVar].to_numpy()

        met = prop.shape[0] * ['']
        if metVar != '':
            met = tab[metVar].to_numpy()

        # [pre, tem, prop, fid, met]
        data = np.hstack((np.vstack(pre), np.vstack(tem), np.vstack(prop), np.vstack(fid), np.vstack(met)))

        return data

    def checkData(self, data, defaultPressure):
        # check pressure
        data = checkPressure(self.preCol, data, defaultPressure)

        # check temperature
        data = checkTemperature(self.temCol, data)

        # check property value
        propVar = self.var
        checkProperty(self.propCol, data, propVar)

        return data

    def convertTemperature(self, col, data):
        tem_convert = self.tem_convert
        if tem_convert != 0.0:
            data[:, col] = tem_convert + data[:, col].astype(np.float)

    def convertProperty(self, col, data):
        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, col] = prop_convert * data[:, col]


class dbsCompound(DbsRelation):
    def getData(self, tab, defaultPressure):
        raise myExceptions.MethodNotImplemented('Method not implemented: dbsRelation::getData{}')

    def getLabel(self):
        raise myExceptions.MethodNotImplemented('Method not implemented: dbsRelation::getLabel{}')


class dbsDensity(DbsRelation):
    def getLabel(self):
        return r'$\rho_{liq} \/ [kg \cdot m^{-1}]$'


class dbsVaporizationEnthalpy(DbsRelation):
    def getLabel(self):
        return r'$\Delta H_{vap} \/ [kJ \cdot mol^{-1}]$'


class dbsVaporizationEnthalpyAtBoilingPoint(DbsRelation):
    def getLabel(self):
        return r'$\Delta H_{vap} \/ [kJ \cdot mol^{-1}]$'


class dbsMeltingPoint(DbsRelation):
    def convertProperty(self, col, data):
        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, col] = prop_convert + data[:, col].astype(np.float)

    def getLabel(self):
        return r'$Tm [K]$'


class dbsBoilingPoint(DbsRelation):
    def convertProperty(self, col, data):
        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, col] = prop_convert + data[:, col].astype(np.float)

    def getLabel(self):
        return r'$Tb [K]$'


class dbsCriticalTemperature(DbsRelation):
    def convertProperty(self, col, data):
        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, col] = prop_convert + data[:, col]

    def getLabel(self):
        return r'$Tc [K]$'


class dbsVaporPressure(DbsRelation):
    def getLabel(self):
        return r'$v_P [bar]$'


class dbsSurfaceTension(DbsRelation):
    def getLabel(self):
        return r'$\gamma \, [mN \cdot m^{-1}]$'


class dbsIsothermalCompressibility(DbsRelation):
    def getLabel(self):
        return r'$\kappa \, [10^{-5} bar^{-1}]$'


class dbsThermalExpansionCoefficient(DbsRelation):
    def getLabel(self):
        return r'$\alpha \, [10^{-4} \, K^{-1}]$'


class dbsHeatCapacityAtConstantPressure(DbsRelation):
    def getLabel(self):
        return r'$c_P \, [JK^{-1} \, mol^{-1}]$'


class dbsPermittivity(DbsRelation):
    def getLabel(self):
        return r'$\epsilon$'


class dbsSelfDiffusionCoefficient(DbsRelation):
    def getLabel(self):
        return r'$D \, [10^{-9} \, m^2 \, s^{-1}]$'


class dbsViscosity(DbsRelation):
    def getLabel(self):
        return r'$\eta \, [10^{-3} \, Pa \cdot s]$'
