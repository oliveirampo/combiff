from abc import ABC, abstractmethod
import json
import sys
import os


import dbsRelation
import myExceptions
import dbsIO
import Equation


class abstractFactoryRelation(ABC):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, eqn, fid, met,
                 prop_convert, tem_convert, pre_convert, marker, color):
        self.larsCode = larsCode
        self.prop = prop
        self.rel = rel
        self.var = var
        self.col = col
        #
        self.pre = pre
        self.tem = tem
        #
        self.eqn = eqn
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

    @abstractmethod
    def createRelation(self):
        pass

    @abstractmethod
    def createParser(self):
        pass

    @abstractmethod
    def createEquations(self, tab):
        pass


class createRelationCompound(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsCompound(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)


    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


    def createEquations(self, tab):
        return [Equation.nullEquation(tab)]


class createRelationDensity(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsDensity(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)


    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


    def createEquations(self, tab):
        if self.eqn == '':
            return [Equation.nullEquation(tab)]
        elif self.eqn == 'dns_1':
            return [Equation.dnsEquation1(tab)]
        else:
            typ = "\'\'".format(self.eqn)
            raise myExceptions.EquationNotImplemented(typ)


class createRelationVaporizationEnthalpy(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsVaporizationEnthalpy(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)


    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


    def createEquations(self, tab):
        if self.eqn == '':
            return [Equation.nullEquation(tab)]
        elif self.eqn == 'hvp_1':
            return [Equation.hvpEquation1(tab)]
        else:
            raise myExceptions.EquationNotImplemented(type)


class createRelationVaporizationEnthalpyAtBoilingPoint(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsVaporizationEnthalpyAtBoilingPoint(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)


    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


    def createEquations(self, tab):
        return [Equation.nullEquation(tab)]


class createRelationBoilingPoint(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsBoilingPoint(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)


    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


    def createEquations(self, tab):
        return [Equation.nullEquation(tab)]


class createRelationMeltingPoint(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsMeltingPoint(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)


    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


    def createEquations(self, tab):
        return [Equation.nullEquation(tab)]


class createRelationCriticalTemperature(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsCriticalTemperature(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)


    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


    def createEquations(self, tab):
        return [Equation.nullEquation(tab)]


class createRelationVaporPressure(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsVaporPressure(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)


    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


    def createEquations(self, tab):
        return [Equation.pvpEquation1(tab)]


class createRelationSurfaceTension(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsSurfaceTension(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)

    def createEquations(self, tab):
        if self.eqn == '':
            return [Equation.nullEquation(tab)]
        elif self.eqn == 'gam_1':
            return [Equation.gamEquation1(tab)]
        else:
            typ = "\'\'".format(self.eqn)
            raise myExceptions.EquationNotImplemented(typ)


class createRelationIsothermalCompressibility(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsIsothermalCompressibility(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)

    def createEquations(self, tab):
        return [Equation.nullEquation(tab)]


class createRelationThermalExpansionCoefficient(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsThermalExpansionCoefficient(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)

    def createEquations(self, tab):
        if self.eqn == '':
            return [Equation.nullEquation(tab)]
        elif self.eqn == 'alp_1':
            return [Equation.alpEquation1(tab)]
        else:
            typ = "\'\'".format(self.eqn)
            raise myExceptions.EquationNotImplemented(typ)


class createRelationHeatCapacityAtConstantPressure(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsHeatCapacityAtConstantPressure(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)

    def createEquations(self, tab):
        if self.eqn == '':
            return [Equation.nullEquation(tab)]
        elif self.eqn == 'hcp_123':
            equations = []

            for idx, row in tab.iterrows():
                df = row.to_frame().T
                equations.append(Equation.hcpEquation1(df))
                equations.append(Equation.hcpEquation2(df))
                equations.append(Equation.hcpEquation3(df))
            # return [Equation.hcpEquation1(tab), Equation.hcpEquation2(tab), Equation.hcpEquation3(tab)]
            return equations
        else:
            typ = "\'\'".format(self.eqn)
            raise myExceptions.EquationNotImplemented(typ)


class createRelationPermittivity(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsPermittivity(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)

    def createEquations(self, tab):
        return [Equation.nullEquation(tab)]


class createRelationSelfDiffusionCoefficient(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsSelfDiffusionCoefficient(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)

    def createEquations(self, tab):
        return [Equation.nullEquation(tab)]


class createRelationViscosity(abstractFactoryRelation):
    def createRelation(self):
        return dbsRelation.dbsViscosity(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.fid, self.met, self.prop_convert, self.tem_convert, self.pre_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)

    def createEquations(self, tab):
        return [Equation.nullEquation(tab)]


class dbsEntry:
    classes = {'cpd': createRelationCompound, 'dns': createRelationDensity, 'hvp': createRelationVaporizationEnthalpy,
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
        try:
            factory = dbsEntry.classes[prop](larsCode, prop, rel, var, col, pre, tem, eqn, fid, met, prop_convert,
                                             tem_convert, pre_convert, marker, color)
        except TypeError as err:
            raise myExceptions.MethodNotImplemented(err.args[0])
        return factory

    def addFactory(self, factory):
        prop = factory.prop
        if prop in self.factories:
            s = '{}: Relation already added.'.format(prop)
            print(s)
            sys.exit(666)
        self.factories[prop] = factory

    def getFactory(self, rel):
        if rel not in self.factories:
            s = 'dbs.json({}) or relation not implemented'.format(self.larsCode)
            raise myExceptions.NoKey(rel, s)
        return self.factories[rel]


def dbsEntryDecoder(obj):
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
                if 'color' in  obj[prop]:
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
    if not os.path.exists(fileName):
        raise myExceptions.NoFile(fileName)

    with open(fileName) as jsonFile:
        dbsEntries = json.load(jsonFile, object_hook=dbsEntryDecoder)

        return dbsEntries


