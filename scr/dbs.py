from abc import ABC, abstractmethod
import json
import sys
import os


import dbsRelation
import myExceptions
from scr import dbsIO


class abstractFactoryRelation(ABC):
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
    def createRelation(self):
        pass

    @abstractmethod
    def createParser(self):
        pass


class createRelationCompound(abstractFactoryRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(createRelationCompound, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)

    def createRelation(self):
        return dbsRelation.dbsCompound(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.prop_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


class createRelationDensity(abstractFactoryRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(createRelationDensity, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)

    def createRelation(self):
        return dbsRelation.dbsDensity(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.prop_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


class createRelationVaporizationEnthalpy(abstractFactoryRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(createRelationVaporizationEnthalpy, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)

    def createRelation(self):
        return dbsRelation.dbsVaporizationEnthalpy(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.prop_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


class createRelationVaporizationEnthalpyAtBoilingPoint(abstractFactoryRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(createRelationVaporizationEnthalpyAtBoilingPoint, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)

    def createRelation(self):
        return dbsRelation.dbsVaporizationEnthalpyAtBoilingPoint(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.prop_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


class createRelationBoilingPoint(abstractFactoryRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(createRelationBoilingPoint, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)

    def createRelation(self):
        return dbsRelation.dbsBoilingPoint(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.prop_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


class createRelationCriticalTemperature(abstractFactoryRelation):
    def __init__(self, larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        super(createRelationCriticalTemperature, self).__init__(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)

    def createRelation(self):
        return dbsRelation.dbsCriticalTemperature(self.larsCode, self.prop, self.rel, self.var, self.col, self.pre, self.tem, self.prop_convert, self.marker, self.color)

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


class dbsEntry():
    path = '/home/marina/PycharmProjects/combiff/wrkDir/tmp'
    propCod = {'dns': ['YA14.6', 'RU18.1', 'FR06.6'], 'hvp': ['YA14.6', 'AC16.1'], 'hvb': ['YA14.6'], 'blp': ['YA14.6'],
               'tem_cri': ['YA14.6']}
    classes = {'cpd': createRelationCompound, 'dns': createRelationDensity, 'hvp': createRelationVaporizationEnthalpy,
               'hvb': createRelationVaporizationEnthalpyAtBoilingPoint, 'blp': createRelationBoilingPoint,
               'tem_cri': createRelationCriticalTemperature}

    def __init__(self, larsCode):
        self.larsCode = larsCode
        self.factories = {}

    def __str__(self):
        s = '\tLARSCODE: {}\n\tRELATIONS:'.format(self.larsCode)
        for rel in self.factories:
            s = '{} {},'.format(s, rel)
        s = s[:-1]
        return s

    def createFactory(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color):
        factory = dbsEntry.classes[prop](larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)
        return factory

    def addFactory(self, factory):
        prop = factory.prop
        if prop in self.factories:
            s = '{}: Relation already added.'.format(prop)
            print(s)
            sys.exit(666)
        self.factories[prop] = factory

    def getFactory(self, rel):
        return self.factories[rel]


def dbsEntryDecoder(obj):
    if '__type__' in obj and obj['__type__'] == 'dbsEntry':
        larsCode = obj['larsCode']
        entry = dbsEntry(larsCode)

        classes = dbsEntry.classes
        for prop in classes:
            if prop in obj:
                rel = obj[prop]['rel']
                var = obj[prop]['var']
                col = obj[prop]['col']

                pre = ''
                if 'pre' in obj[prop]:
                    pre = obj[prop]['pre']

                tem = ''
                if 'tem' in obj[prop]:
                    tem = obj[prop]['tem']

                prop_convert = 1.0
                if 'prop_convert' in obj[prop]:
                    prop_convert = obj[prop]['prop_convert']
                    prop_convert = int(prop_convert)

                marker = ''
                if 'marker' in  obj[prop]:
                    marker = obj[prop]['marker']

                color = ''
                if 'color' in  obj[prop]:
                    color = obj[prop]['color']

                factory = dbsEntry.createFactory(larsCode, prop, rel, var, col, pre, tem, prop_convert, marker, color)
                entry.addFactory(factory)

        return entry
    else:
        return obj


def getDbsEntries(fileName):
    if not os.path.exists(fileName):
        raise myExceptions.NoFile(fileName)

    with open(fileName) as jsonFile:
        dbsEntries = json.load(jsonFile, object_hook=dbsEntryDecoder)
        # for entry in dbsEntries: print(dbsEntries[entry])

        return dbsEntries


