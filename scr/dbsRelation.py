from abc import ABC, abstractmethod
import numpy as np
import sys


import myExceptions


class DbsRelation(ABC):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color):
        self.larsCode = larsCode
        self.prop = prop
        self.rel = rel
        self.var = var
        self.col = col
        #
        self.pre = pre
        self.tem = tem
        #
        self.prop_convert = prop_convert
        self.tem_convert = tem_convert
        self.pre_convert = pre_convert
        #
        self.marker = marker
        self.color = color


    @abstractmethod
    def getData(self, tab, defaultPressure):
        pass


    def getDataHelper(self, tab, defaultPressure):
        columns = self.col.strip().split()

        if columns[1] != 'rec':
            print('Problem with col sequence for {}_{}'.format(self.rel, self.larsCode))
            sys.exit(666)

        preVar = self.pre
        temVar = self.tem
        propVar = self.var
        # print(preVar, temVar, propVar)

        if (preVar != '') and (temVar != ''):
            data = tab[[preVar, temVar, propVar]]
            data = data.to_numpy()
        elif (preVar == '') and (temVar == ''):
            data = tab[[propVar]]
            data = data.to_numpy()
        else:
            data = tab[[temVar, propVar]]
            data = data.to_numpy()

            pre = data.shape[0] * [[defaultPressure]]
            data = np.hstack((pre, data))

        return data


    # if pre == '%' -> pre = 1.0 kPa
    def checkPressure(self, col, data, defaultPressure):
        for i in range(data[:, col].shape[0]):
            try:
                data[i, col] = float(data[i, col])

            except ValueError:
                val = data[i, col]
                if val == '%':
                    data[i, col] = defaultPressure
                else:
                    raise myExceptions.WrongProperty('pressure', val)


    # if tem == '%' -> remove data
    def checkTemperature(self, col, data):
        for i in range(data[:, col].shape[0]):
            try:
                data[i, col] = float(data[i, col])

            except ValueError:
                val = data[i, col]
                raise myExceptions.WrongProperty('temperature', val)


    # if val == '%' -> remove data
    def checkProperty(self, col, data, propVar):
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
                    raise myExceptions.WrongProperty(propVar, val)


    def checkTransitionTemperature(self, data):
        # check property value
        col = 0
        propVar = self.var
        self.checkProperty(col, data, propVar)

        # convert temperature
        tem_convert = self.tem_convert
        if tem_convert != 0:
            data = data + tem_convert


class dbsCompound(DbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color):
        super(dbsCompound, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color)


    def getData(self, tab, defaultPressure):
        raise myExceptions.MethodNotImplemented('Method not implemented: getData')


class dbsDensity(DbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color):
        if prop_convert == '':
            prop_convert = 1.0
        super(dbsDensity, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color)


    def getData(self, tab, defaultPressure):
        data = super(dbsDensity, self).getDataHelper(tab, defaultPressure)

        # check pressure
        col = 0
        self.checkPressure(col, data, defaultPressure)

        # check temperature
        col = 1
        self.checkTemperature(col, data)

        # check property value
        col = 2
        propVar = self.var
        self.checkProperty(col, data, propVar)

        # convert temperature
        tem_convert = self.tem_convert
        if tem_convert != 0.0:
            data[:, 1] = tem_convert + data[:, 1]

        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, 2] = prop_convert * data[:, 2]

        # print(data)
        return data


class dbsVaporizationEnthalpy(DbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color):
        if prop_convert == '':
            prop_convert = 1.0
        super(dbsVaporizationEnthalpy, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color)


    def getData(self, tab, defaultPressure):
        data = super(dbsVaporizationEnthalpy, self).getDataHelper(tab, defaultPressure)

        # check pressure
        col = 0
        self.checkPressure(col, data, defaultPressure)

        # check temperature
        col = 1
        self.checkTemperature(col, data)

        # check property value
        col = 2
        propVar = self.var
        self.checkProperty(col, data, propVar)

        # convert temperature
        tem_convert = self.tem_convert
        if tem_convert != 0.0:
            data[:, 1] = tem_convert + data[:, 1]

        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, 2] = prop_convert * data[:, 2]

        # print(data)
        return data


class dbsVaporizationEnthalpyAtBoilingPoint(DbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color):
        if prop_convert == '':
            prop_convert = 1.0
        super(dbsVaporizationEnthalpyAtBoilingPoint, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color)


    def getData(self, tab, defaultPressure):
        data = super(dbsVaporizationEnthalpyAtBoilingPoint, self).getDataHelper(tab, defaultPressure)

        # check pressure
        col = 0
        self.checkPressure(col, data, defaultPressure)

        # check temperature
        col = 1
        self.checkTemperature(col, data)

        # check property value
        col = 2
        propVar = self.var
        self.checkProperty(col, data, propVar)

        # convert temperature
        tem_convert = self.tem_convert
        if tem_convert != 0.0:
            data[:, 1] = tem_convert + data[:, 1]

        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, 2] = prop_convert * data[:, 2]

        return data


class dbsMeltingPoint(DbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color):
        if prop_convert == '':
            prop_convert = 0.0
        super(dbsMeltingPoint, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color)


    def getData(self, tab, defaultPressure):
        data = super(dbsMeltingPoint, self).getDataHelper(tab, defaultPressure)
        super(dbsMeltingPoint, self).checkTransitionTemperature(data)
        return data


class dbsBoilingPoint(DbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color):
        if prop_convert == '':
            prop_convert = 0.0
        super(dbsBoilingPoint, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color)


    def getData(self, tab, defaultPressure):
        data = super(dbsBoilingPoint, self).getDataHelper(tab, defaultPressure)
        super(dbsBoilingPoint, self).checkTransitionTemperature(data)
        return data


class dbsCriticalTemperature(DbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color):
        if prop_convert == '':
            prop_convert = 0.0
        super(dbsCriticalTemperature, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, tem_convert, pre_convert, marker, color)


    def getData(self, tab, defaultPressure):
        data = super(dbsCriticalTemperature, self).getDataHelper(tab, defaultPressure)
        super(dbsCriticalTemperature, self).checkTransitionTemperature(data)
        return data




