"""Module that handles DBS relations.

Classes:
    DbsRelation
    dbsCompound(DbsRelation)
    dbsDensity(DbsRelation)
    dbsVaporizationEnthalpy(DbsRelation)
    dbsVaporizationEnthalpyAtBoilingPoint(DbsRelation)
    dbsMeltingPoint(DbsRelation)
    dbsBoilingPoint(DbsRelation)
    dbsCriticalTemperature(DbsRelation)
    dbsVaporPressure(DbsRelation)
    dbsSurfaceTension(DbsRelation)
    dbsIsothermalCompressibility(DbsRelation)
    dbsThermalExpansionCoefficient(DbsRelation)
    dbsHeatCapacityAtConstantPressure(DbsRelation)
    dbsPermittivity(DbsRelation)
    dbsSelfDiffusionCoefficient(DbsRelation)
    dbsViscosity(DbsRelation)

Methods:
    checkPressure(col, data, defaultPressure)
    checkTemperature(col, data)
    checkProperty(col, data, propVar)
    getUnit(prop)
"""

import numpy as np
import sys


from scr.base import myExceptions


def checkPressure(col, data, defaultPressure):
    """Checks if pressure values are consistent.
    Of pre == '%' -> pre = 1.0 kPa

    :param col: (str) Name of columns.
    :param data: (pandas DataFrame) Table with data.
    :param defaultPressure: (float) Default pressure.
    :return:
        data: (pandas DataFrame) Table with data.
    """

    rm = []

    for i in range(data[:, col].shape[0]):
        try:
            data[i, col] = float(data[i, col])

        except ValueError:
            val = data[i, col]
            if val == '%':
                data[i, col] = defaultPressure

            elif val == 'sat._liquid':
                rm.append(i)

            elif val == 'sat__liq':
                rm.append(i)

            else:
                raise myExceptions.WrongProperty('pressure', val, 'dbsRelation::checkPressure()')

    for idx in reversed(rm):
        data = np.delete(data, idx, 0)

    return data


def checkTemperature(col, data):
    """Checks if temperature values are consistent.
    If tem == '%' -> remove row.

    :param col: (str) Name of columns.
    :param data: (pandas DataFrame) Table with data.
    :return:
        data: (pandas DataFrame) Table with data.
    """

    rm = []

    for i in range(data[:, col].shape[0]):
        try:
            data[i, col] = float(data[i, col])

        except ValueError:
            val = data[i, col]
            if val == '%':
                rm.append(i)
            elif val == 'sat._liquid':
                rm.append(i)
            elif val.endswith('?'):
                data[i, col] = val[:-1]
            else:
                print(data[i, col])
                raise myExceptions.WrongProperty('temperature', val, 'dbsRelation:checkTemperature()')

    for idx in reversed(rm):
        data = np.delete(data, idx, 0)

    return data


def checkProperty(col, data, propVar):
    """Checks if property values are consistent.
    If val == '%' -> remove data.

    :param col: (str) Name of columns.
    :param data: (pandas DataFrame) Table with data.
    :param propVar: (str) Code of relation.
    :return:
        data: (pandas DataFrame) Table with data.
    """

    rm = []

    for i in range(data[:, col].shape[0]):
        try:
            data[i, col] = float(data[i, col])

        except ValueError:
            val = data[i, col]
            if '+/-' in val:
                val = val.split('+/-')
                val = float(val[0])
                data[i, col] = val

            elif '<' in val:
                rm.append(i)

            elif ',' in val:
                val = val.replace(',', '.')
                val = float(val[0])
                data[i, col] = val

            else:
                print(data[i, col])
                raise myExceptions.WrongProperty(propVar, val, 'dbsRelation::checkProperty()')

    for idx in reversed(rm):
        data = np.delete(data, idx, 0)

    return data


def getUnit(prop):
    """Returns unit property.

    :param prop: (str) Property code.
    :return: (str)
    """

    if prop == 'dns':
        return dbsDensity.getUnit()
    elif prop == 'hvp':
        return dbsVaporizationEnthalpy.getUnit()
    elif prop == 'gam':
        return dbsSurfaceTension.getUnit()
    elif prop == 'kap':
        return dbsIsothermalCompressibility.getUnit()
    elif prop == 'alp':
        return dbsThermalExpansionCoefficient.getUnit()
    elif prop == 'hcp':
        return dbsHeatCapacityAtConstantPressure.getUnit()
    elif prop == 'eps':
        return dbsPermittivity.getUnit()
    elif prop == 'diffus':
        return dbsSelfDiffusionCoefficient.getUnit()
    elif prop == 'etd':
        return dbsViscosity.getUnit()
    else:
        print('No such property: {}'.format(prop))
        sys.exit(1)


class DbsRelation:
    """DbsRelation object.

    Attributes:
        larsCode: (str) LARS code.
        prop: (prop) Property code.
        rel: (str) Relation code.
        var: (str) Code for property in relation table.
        col: (str) Name of columns.
        pre: (str) Code for pressure in relation table.
        tem: (str) Code for temperature in relation table.
        fid: (str) Code for annotation in relation table.
        met: (str) Code for method in relation table.
        prop_convert: (str) Constant to convert property unit.
        tem_convert: (str) Constant to convert temperature unit.
        pre_convert: (str) Constant to convert pressure unit.
        marker: (str) Marker to be used in plot.
        color: (str) Color to be used in plot.

        preCol: (int, 0)
        temCol: (int, 1)
        propCol: (int, 2)
        fidCol: (int, 3)
        metCol: (int, 4)

    Functions:
        getUnit()
        getMarker(self)
        getColor(self)
        getValues(data)
        getData(tab, defaultPressure)
        getDataHelper(tab, defaultPressure)
        checkData(data, defaultPressure)
        convertTemperature(col, data)
        convertProperty(col, data)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):
        """Constructs all the necessary attributes for this object.

        :param larsCode: (str) LARS code.
        :param prop: (prop) Property code.
        :param rel: (str) Relation code.
        :param var: (str) Code for property in relation table.
        :param col: (str) Name of columns.
        :param pre: (str) Code for pressure in relation table.
        :param tem: (str) Code for temperature in relation table.
        :param fid: (str) Code for annotation in relation table.
        :param met: (str) Code for method in relation table.
        :param prop_convert: (str) Constant to convert property unit.
        :param tem_convert: (str) Constant to convert temperature unit.
        :param pre_convert: (str) Constant to convert pressure unit.
        :param marker: (str) Marker to be used in plot.
        :param color: (str) Color to be used in plot.
        """

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

    @staticmethod
    def getUnit():
        """Returns unit of property in relation."""

        pass

    def getMarker(self):
        """Returns marker to be used in plot."""

        return self.marker

    def getColor(self):
        """Returns color to be used in plot."""

        return self.color

    def getValues(self, data):
        """Returns values from table.

        :param data: (pandas DataFrame) Table with data.
        :return:
            pre: (numpy array) Pressure.
            tem: (numpy array) Temperature.
            val: (numpy array) Property value.
            fid: (numpy array) Annotation.
            met: (numpy array) Method.
        """

        pre = data[:, self.preCol]
        tem = data[:, self.temCol]
        val = data[:, self.propCol]
        fid = data[:, self.fidCol]
        met = data[:, self.metCol]

        pre = pre.astype(np.float)
        tem = tem.astype(np.float)
        val = val.astype(np.float)

        return pre, tem, val, fid, met

    def getData(self, tab, defaultPressure):
        """Returns data (pressure, temperature and value) for the current relation.

        :param tab: (pandas DataFrame) Table with data.
        :param defaultPressure: (float) Default pressure.
        :return:
            data: (pandas DataFrame) Table with data.
        """

        data = self.getDataHelper(tab, defaultPressure)
        data = self.checkData(data, defaultPressure)

        col = self.temCol
        self.convertTemperature(col, data)

        col = self.propCol
        self.convertProperty(col, data)

        return data

    def getDataHelper(self, tab, defaultPressure):
        """Helper function to extract data from table.

        :param tab: (pandas DataFrame) Table with all data.
        :param defaultPressure: (float) Default pressure.
        :return:
            data: (pandas DataFrame) Table with data [pressure, temperature, property value, annotation, method]
        """

        columns = self.col.strip().split()

        if columns[1] != 'rec':
            print('Problem with col sequence for {}_{}'.format(self.rel, self.larsCode))
            sys.exit(666)

        preVar = self.pre
        temVar = self.tem
        propVar = self.var
        fidVar = self.fid
        metVar = self.met

        if propVar not in tab:
            print(tab.columns)
            print(self.larsCode)
            sys.exit(123)
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
        """Checks if data is consistent.

        :param data: (numpy.ndarray) Data.
        :param defaultPressure: (float) Default pressure.
        :return:
            (numpy.ndarray) Data.
        """

        # check pressure
        data = checkPressure(self.preCol, data, defaultPressure)

        # check temperature
        data = checkTemperature(self.temCol, data)

        # check property value
        propVar = self.var
        data = checkProperty(self.propCol, data, propVar)

        return data

    def convertTemperature(self, col, data):
        """Converts temperature unit.

        :param col: (str) Column names.
        :param data: (numpy.ndarray) Data.
        """

        tem_convert = self.tem_convert
        if tem_convert != 0.0:
            data[:, col] = tem_convert + data[:, col].astype(np.float)

    def convertProperty(self, col, data):
        """Converts pressure unit.

        :param col: (str) Column names.
        :param data: (numpy.ndarray) Data.
        :return:
        """

        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, col] = prop_convert * data[:, col].astype(np.float)


class dbsCompound(DbsRelation):
    """Implements DBS compound relation for a given source (LARS code)."""

    def getData(self, tab, defaultPressure):
        raise myExceptions.MethodNotImplemented('Method not implemented: dbsRelation::getData{}')

    @staticmethod
    def getUnit():
        raise myExceptions.MethodNotImplemented('Method not implemented: dbsRelation::getUnit{}')


class dbsDensity(DbsRelation):
    """Implements DBS density relation for a given source (LARS code)."""

    @staticmethod
    def getUnit():
        return r'$\rho_{liq} \/ [kg \cdot m^{-1}]$'


class dbsVaporizationEnthalpy(DbsRelation):
    """Implements DBS vaporization enthalpy relation for a given source (LARS code)."""

    @staticmethod
    def getUnit():
        return r'$\Delta H_{vap} \/ [kJ \cdot mol^{-1}]$'


class dbsVaporizationEnthalpyAtBoilingPoint(DbsRelation):
    """Implements DBS vaporization enthalpy at constant pressure relation for a given source (LARS code)."""

    @staticmethod
    def getUnit():
        return r'$\Delta H_{vap} \/ [kJ \cdot mol^{-1}]$'


class dbsMeltingPoint(DbsRelation):
    """Implements DBS melting point relation for a given source (LARS code)."""

    def convertProperty(self, col, data):
        """Converts temperature unit.

        :param col: (str) Column names.
        :param data: (numpy.ndarray) Data.
        """

        if data.shape[0] == 0:
            return

        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, col] = prop_convert + data[:, col].astype(np.float)

    @staticmethod
    def getUnit():
        return r'$Tm [K]$'


class dbsBoilingPoint(DbsRelation):
    """Implements DBS boiling point relation for a given source (LARS code)."""

    def convertProperty(self, col, data):
        """Converts temperature unit.

        :param col: (str) Column names.
        :param data: (numpy.ndarray) Data.
        """

        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, col] = prop_convert + data[:, col].astype(np.float)

    @staticmethod
    def getUnit():
        return r'$Tb [K]$'


class dbsCriticalTemperature(DbsRelation):
    """Implements DBS critical temperature relation for a given source (LARS code)."""

    def convertProperty(self, col, data):
        """Converts temperature unit.

        :param col: (str) Column names.
        :param data: (numpy.ndarray) Data.
        """

        prop_convert = self.prop_convert
        if prop_convert != 1:
            data[:, col] = prop_convert + data[:, col]

    @staticmethod
    def getUnit():
        return r'$Tc [K]$'


class dbsVaporPressure(DbsRelation):
    """Implements DBS vapor pressure relation for a given source (LARS code)."""

    @staticmethod
    def getUnit():
        return r'$v_P [bar]$'


class dbsSurfaceTension(DbsRelation):
    """Implements DBS surface tension relation for a given source (LARS code)."""
    @staticmethod
    def getUnit():
        return r'$\gamma \, [mN \cdot m^{-1}]$'


class dbsIsothermalCompressibility(DbsRelation):
    """Implements DBS isothermal compressibility relation for a given source (LARS code)."""

    @staticmethod
    def getUnit():
        # return r'$\kappa \, [10^{-5} bar^{-1}]$'
        return r'$\kappa \, [bar^{-1}]$'


class dbsThermalExpansionCoefficient(DbsRelation):
    """Implements DBS thermal expansion coefficient relation for a given source (LARS code)."""

    @staticmethod
    def getUnit():
        return r'$\alpha \, [10^{-4} \, K^{-1}]$'


class dbsHeatCapacityAtConstantPressure(DbsRelation):
    """Implements DBS heat capacity at constant pressure relation for a given source (LARS code)."""

    @staticmethod
    def getUnit():
        return r'$c_P \, [JK^{-1} \, mol^{-1}]$'


class dbsPermittivity(DbsRelation):
    """Implements DBS permittivity relation for a given source (LARS code)."""

    @staticmethod
    def getUnit():
        return r'$\epsilon$'


class dbsSelfDiffusionCoefficient(DbsRelation):
    """Implements DBS self diffusion coefficient relation for a given source (LARS code)."""

    @staticmethod
    def getUnit():
        return r'$D \, [10^{-9} \, m^2 \, s^{-1}]$'


class dbsViscosity(DbsRelation):
    """Implements DBS viscosity relation for a given source (LARS code)."""

    @staticmethod
    def getUnit():
        return r'$\eta \, [10^{-3} \, Pa \cdot s]$'
