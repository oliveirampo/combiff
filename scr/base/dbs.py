from abc import ABC, abstractmethod
import pandas as pd
import json
import sys
import os


from scr.base import dbsRelation
from scr.base import myExceptions
from scr.base import dbsIO
from scr.base import equation


class abstractFactoryRelation(ABC):
    """Abstract class to factory to create relations.

    Attributes:
        _prop: (str) Property code.
        _eqn: (str) Code for equation in relation table.
        _var: (str) Column name in table.
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, prop, var, eqn):
        """Constructs all the necessary attributes for this object.

        :param prop: (str) Property code.
        :param var: (str) Column name in table.
        :param eqn: (str) Code for equation in relation table.
        """

        self._prop = prop
        self._var = var
        self._eqn = eqn

    @abstractmethod
    def getRelation(self):
        """Returns created DBS relation."""
        pass

    @abstractmethod
    def getParser(self):
        """Returns created relationParser"""
        pass

    @abstractmethod
    def createEquations(self, tab):
        """Creates and returns list of equations."""
        pass

    @property
    def prop(self):
        return self._prop

    @property
    def eqn(self):
        return self._eqn

    @property
    def var(self):
        return self._var


class createRelationCompound(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsCompound(larsCode, prop, rel, var, col, pre, tem, fid, met, prop_convert,
                                                    tem_convert, pre_convert, marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        return [equation.nullEquation(tab)]


class createRelationDensity(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsDensity(larsCode, prop, rel, var, col, pre, tem, fid, met, prop_convert,
                                                   tem_convert, pre_convert, marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        if self.eqn == '':
            return [equation.nullEquation(tab)]
        elif self.eqn == 'dns_1':
            equations = []

            for i in range(tab.shape[0]):
                newTab = tab.iloc[i].to_frame().T
                equations.append(equation.dnsEquation1(newTab))

            return equations

        else:
            typ = "\'\'".format(self.eqn)
            raise myExceptions.EquationNotImplemented(typ)


class createRelationVaporizationEnthalpy(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsVaporizationEnthalpy(larsCode, prop, rel, var, col, pre, tem, fid, met,
                                                                prop_convert, tem_convert, pre_convert, marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        if self.eqn == '':
            return [equation.nullEquation(tab)]
        elif self.eqn == 'hvp_1':
            equations = []

            for i in range(tab.shape[0]):
                newTab = tab.iloc[i].to_frame().T
                equations.append(equation.hvpEquation1(newTab))

            return equations
        else:
            raise myExceptions.EquationNotImplemented(type)


class createRelationVaporizationEnthalpyAtBoilingPoint(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsVaporizationEnthalpyAtBoilingPoint(larsCode, prop, rel, var, col, pre, tem,
                                                                              fid, met, prop_convert, tem_convert,
                                                                              pre_convert, marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        return [equation.nullEquation(tab)]


class createRelationBoilingPoint(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsBoilingPoint(larsCode, prop, rel, var, col, pre, tem, fid, met, prop_convert,
                                                        tem_convert, pre_convert, marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        return [equation.nullEquation(tab)]


class createRelationMeltingPoint(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):
        self._dbsRelation = dbsRelation.dbsMeltingPoint(larsCode, prop, rel, var, col, pre, tem, fid, met, prop_convert,
                                                        tem_convert, pre_convert, marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        return [equation.nullEquation(tab)]


class createRelationCriticalTemperature(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsCriticalTemperature(larsCode, prop, rel, var, col, pre, tem, fid, met,
                                                               prop_convert, tem_convert, pre_convert, marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        return [equation.nullEquation(tab)]


class createRelationVaporPressure(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsVaporPressure(larsCode, prop, rel, var, col, pre, tem, fid, met,
                                                         prop_convert, tem_convert, pre_convert, marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        equations = []

        for i in range(tab.shape[0]):
            newTab = tab.iloc[i].to_frame().T
            equations.append(equation.pvpEquation1(newTab))

        return equations


class createRelationSurfaceTension(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsSurfaceTension(larsCode, prop, rel, var, col, pre, tem, fid, met,
                                                          prop_convert, tem_convert, pre_convert, marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        if self.eqn == '':
            return [equation.nullEquation(tab)]
        elif self.eqn == 'gam_1':
            return [equation.gamEquation1(tab)]
        else:
            typ = "\'\'".format(self.eqn)
            raise myExceptions.EquationNotImplemented(typ)


class createRelationIsothermalCompressibility(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsIsothermalCompressibility(larsCode, prop, rel, var, col, pre, tem, fid, met,
                                                                     prop_convert, tem_convert, pre_convert, marker,
                                                                     color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        return [equation.nullEquation(tab)]


class createRelationThermalExpansionCoefficient(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsThermalExpansionCoefficient(larsCode, prop, rel, var, col, pre, tem, fid,
                                                                       met, prop_convert, tem_convert, pre_convert,
                                                                       marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        if self.eqn == '':
            return [equation.nullEquation(tab)]
        elif self.eqn == 'alp_1':
            return [equation.alpEquation1(tab)]
        else:
            typ = "\'\'".format(self.eqn)
            raise myExceptions.EquationNotImplemented(typ)


class createRelationHeatCapacityAtConstantPressure(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsHeatCapacityAtConstantPressure(larsCode, prop, rel, var, col, pre, tem, fid,
                                                                          met, prop_convert, tem_convert, pre_convert,
                                                                          marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        if self.eqn == '':
            return [equation.nullEquation(tab)]
        elif self.eqn == 'hcp_123':
            equations = []

            for idx, row in tab.iterrows():
                df = row.to_frame().T
                equations.append(equation.hcpEquation1(df))
                equations.append(equation.hcpEquation2(df))
                equations.append(equation.hcpEquation3(df))

            return equations
        else:
            typ = "\'\'".format(self.eqn)
            raise myExceptions.EquationNotImplemented(typ)


class createRelationPermittivity(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsPermittivity(larsCode, prop, rel, var, col, pre, tem, fid, met,
                                                        prop_convert, tem_convert, pre_convert, marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        if self.eqn == '':
            return [equation.nullEquation(tab)]
        elif self.eqn == 'eps_1':
            a = tab['a_eps'].values[0]
            if a == '%':
                return [equation.nullEquation(tab)]

            return [equation.epsEquation1(tab)]
        else:
            typ = "\'\'".format(self.eqn)
            raise myExceptions.EquationNotImplemented(typ)


class createRelationSelfDiffusionCoefficient(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsSelfDiffusionCoefficient(larsCode, prop, rel, var, col, pre, tem, fid, met,
                                                                    prop_convert, tem_convert, pre_convert, marker,
                                                                    color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        return [equation.nullEquation(tab)]


class createRelationViscosity(abstractFactoryRelation):
    """Implements factory to create specific relation.

    Attributes:
        _dbsRelation: (dbsRelation)
        _parser: (relationParser)
    Functions:
        createRelation()
        createParser()
        createEquations(tab)
    """

    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert, tem_convert, pre_convert,
                 marker, color):

        self._dbsRelation = dbsRelation.dbsViscosity(larsCode, prop, rel, var, col, pre, tem, fid, met,
                                                     prop_convert, tem_convert, pre_convert, marker, color)

        self._parser = dbsIO.parserDefault(larsCode, rel, col)

        super().__init__(prop, var, eqn)

    def getRelation(self):
        """Returns created DBS relation."""
        return self._dbsRelation

    def getParser(self):
        """Returns created relationParser"""
        return self._parser

    def createEquations(self, tab):
        """Creates and returns list of equations."""
        return [equation.nullEquation(tab)]


class dbsEntry:
    """A dbsEntry corresponds to a LARS code and the relations/tables associated with it.
    For each relation/table there will be one factory method to construct the appropriate object.

    Attributes:
        larsCode: (str) LARS code.
        factories: (dict) Dictionary of factories.

    """

    classes = {'cpd': createRelationCompound,
               'dns': createRelationDensity, 'hvp': createRelationVaporizationEnthalpy,
               'hvb': createRelationVaporizationEnthalpyAtBoilingPoint,
               'mlp': createRelationMeltingPoint, 'blp': createRelationBoilingPoint,
               'tem_cri': createRelationCriticalTemperature,
               'pvp': createRelationVaporPressure,
               'gam': createRelationSurfaceTension,
               'kap': createRelationIsothermalCompressibility,
               'alp': createRelationThermalExpansionCoefficient,
               'hcp': createRelationHeatCapacityAtConstantPressure,
               'eps': createRelationPermittivity,
               'diffus': createRelationSelfDiffusionCoefficient,
               'etd': createRelationViscosity}

    def __init__(self, larsCode):
        """Constructs all the necessary attributes for this object.

        :param larsCode: (str) LARS code.
        """

        self.larsCode = larsCode
        self.factories = {}

    def __str__(self):
        s = '\tLARSCODE: {}\n\tRELATIONS:'.format(self.larsCode)
        for rel in self.factories:
            s = '{} {},'.format(s, rel)
        s = s[:-1]
        return s

    @staticmethod
    def createFactory(larsCode, prop, rel, var, col, pre, tem, eqn, fid, met,
                      prop_convert, tem_convert, pre_convert, marker, color):
        """Creates an object that implements abstractFactoryRelation.

        :param larsCode: (str) LARS code.
        :param prop: (prop) Property code.
        :param rel: (str) Relation code.
        :param var: (str) Code for property in relation table.
        :param col: (str) Name of columns.
        :param pre: (str) Code for pressure in relation table.
        :param tem: (str) Code for temperature in relation table.
        :param eqn: (str) Code for equation in relation table.
        :param fid: (str) Code for annotation in relation table.
        :param met: (str) Code for method in relation table.
        :param prop_convert: (str) Constant to convert property unit.
        :param tem_convert: (str) Constant to convert temperature unit.
        :param pre_convert: (str) Constant to convert pressure unit.
        :param marker: (str) Marker to be used in plot.
        :param color: (str) Color to be used in plot.
        :return:
            factory: Implementation of abstractFactoryRelation.
        """

        try:
            factory = dbsEntry.classes[prop](larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert,
                                             tem_convert, pre_convert, marker, color)
        except TypeError as err:
            raise myExceptions.MethodNotImplemented(err.args[0])
        return factory

    def addFactory(self, factory):
        """Adds a factory (object that implements abstractFactoryRelation) to the dictionary of factories.

        :param factory: (object that implements abstractFactoryRelation)
        """

        prop = factory.prop
        if prop in self.factories:
            s = '{}: Relation already added.'.format(prop)
            print(s)
            sys.exit(666)
        self.factories[prop] = factory

    def getFactory(self, rel):
        """Returns factory (object that implements abstractFactoryRelation) given a relation code.

        :param rel: (str) Relation code.
        :return: (object that implements abstractFactoryRelation)
        """

        if rel not in self.factories:
            s = 'dbs.json({}) or relation not implemented'.format(self.larsCode)
            raise myExceptions.NoKey(rel, s)
        return self.factories[rel]


def dbsEntryDecoder(obj):
    """Method to decode a dbsEntry object.

    :param obj:
    :return:
        mol: (dbsEntry)
    """

    if '__type__' in obj and obj['__type__'] == 'dbsEntry':
        larsCode = obj['larsCode']
        entry = dbsEntry(larsCode)

        classes = dbsEntry.classes
        for prop in classes:
            if prop in obj:
                try:
                    rel = obj[prop]['rel']
                    var = obj[prop]['var']
                    col = obj[prop]['col']
                except KeyError as err:
                    raise myExceptions.NoKey(err.args[0], larsCode + "/" + prop)

                pre = ''
                if 'pre' in obj[prop]:
                    pre = obj[prop]['pre']

                tem = ''
                if 'tem' in obj[prop]:
                    tem = obj[prop]['tem']

                prop_convert = 1.0
                if 'prop_convert' in obj[prop]:
                    prop_convert = obj[prop]['prop_convert']
                    prop_convert = float(prop_convert)

                tem_convert = 0.0
                if 'tem_convert' in obj[prop]:
                    tem_convert = obj[prop]['tem_convert']
                    tem_convert = float(tem_convert)

                pre_convert = 1.0
                if 'pre_convert' in obj[prop]:
                    pre_convert = obj[prop]['pre_convert']
                    pre_convert = float(pre_convert)

                eqn = ''
                if 'eqn' in obj[prop]:
                    eqn = obj[prop]['eqn']

                marker = ''
                if 'marker' in obj[prop]:
                    marker = obj[prop]['marker']

                color = ''
                if 'color' in obj[prop]:
                    color = obj[prop]['color']

                fid = ''
                if 'fid' in obj[prop]:
                    fid = obj[prop]['fid']

                met = ''
                if 'method' in obj[prop]:
                    met = obj[prop]['method']

                factory = dbsEntry.createFactory(larsCode, prop, rel, var, col, pre, tem, eqn, fid, met,
                                                 prop_convert, tem_convert, pre_convert, marker, color)
                entry.addFactory(factory)

        return entry
    else:
        return obj


def getDbsEntries(fileName):
    """Reads dbs entries from file.

    :param fileName: (str) File name.
    :return:
        dbsEntries: (dict) Dictionary of DBS entries.
    """

    if not os.path.exists(fileName):
        raise myExceptions.NoFile(fileName)

    with open(fileName) as jsonFile:
        dbsEntries = json.load(jsonFile, object_hook=dbsEntryDecoder)

        return dbsEntries
