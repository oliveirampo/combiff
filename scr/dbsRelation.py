from abc import ABC, abstractmethod
import numpy as np
import sys


class dbsRelation(ABC):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        self.larsCode = larsCode
        self.prop = prop
        self.rel = rel
        self.var = var
        self.col = col
        #
        self.pre = pre
        self.tem = tem
        self.prop_convert = prop_convert
        #
        self.marker = marker
        self.color = color

    @abstractmethod
    def getData(self, tab):
        pass

    def getDataHelper(self, tab):
        columns = self.col.strip().split()

        if columns[1] != 'rec':
            print('Problem with col sequence for {}_{}'.format(self.rel, self.larsCode))
            sys.exit(666)

        preVar = self.pre
        temVar = self.tem
        propVar = self.var
        # print(preVar, temVar, propVar)

        if preVar != '':
            data = tab[[preVar, temVar, propVar]]
            data = data.to_numpy()
        else:
            data = tab[[temVar, propVar]]
            data = data.to_numpy()

            pre = np.ones((data.shape[0], 1))
            data = np.hstack((pre, data))

        return data


class dbsCompound(dbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(dbsCompound, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)


class dbsDensity(dbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(dbsDensity, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)

    def getData(self, tab):
        data = super(dbsDensity, self).getDataHelper(tab)

        # TODO - check for non-float values

        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, 2] = prop_convert * data[:, 2]

        # print(data)
        return data


class dbsVaporizationEnthalpy(dbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(dbsVaporizationEnthalpy, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)

    def getData(self, tab):
        data = super(dbsVaporizationEnthalpy, self).getDataHelper(tab)

        # TODO - check for non-float values
        for i in range(data.shape[0]):
            val = data[i][2]
            try:
                val = float(val)
            except ValueError:
                if '+/-' in val:
                    val = val.split('+/-')[0]
                else:
                    print(val)
                    sys.exit(666)

            data[i][2] = val

        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, 2] = prop_convert * data[:, 2]

        # print(data)
        return data


class dbsVaporizationEnthalpyAtBoilingPoint(dbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(dbsVaporizationEnthalpyAtBoilingPoint, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)


class dbsBoilingPoint(dbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(dbsBoilingPoint, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)


class dbsCriticalTemperature(dbsRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(dbsCriticalTemperature, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)


