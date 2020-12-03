from abc import ABC, abstractmethod
import sys

from scr import dbsIO


class abstractFactoryRelation(ABC):
    def __init__(self, larsCode, prop, rel, var, col):
        self.larsCode = larsCode
        self.prop = prop
        self.rel = rel
        self.var = var
        self.col = col

    @abstractmethod
    def createRelation(self, rel, var, col):
        pass

    @abstractmethod
    def createParser(self):
        pass


class createRelationCompound(abstractFactoryRelation):
    def __init__(self, larsCode, prop, rel, var, col):
        super(createRelationCompound, self).__init__(larsCode, prop, rel, var, col)

    def createRelation(self):
        return dbsCompound()

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


class createRelationDensity(abstractFactoryRelation):
    def __init__(self, larsCode, prop, rel, var, col):
        super(createRelationDensity, self).__init__(larsCode, prop, rel, var, col)

    def createRelation(self):
        return dbsDensity()

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


class createRelationVaporizationEnthalpy(abstractFactoryRelation):
    def __init__(self, larsCode, prop, rel, var, col):
        super(createRelationVaporizationEnthalpy, self).__init__(larsCode, prop, rel, var, col)

    def createRelation(self):
        return dbsVaporizationEnthalpy()

    def createParser(self):
        return dbsIO.parserDefault(self.larsCode, self.rel, self.col)


class dbsRelation(ABC):
    pass


class dbsCompound(dbsRelation):
    pass


class dbsDensity(dbsRelation):
    pass


class dbsVaporizationEnthalpy(dbsRelation):
    pass


class dbsEntry():
    classes = {'cpd': createRelationCompound, 'dns': createRelationDensity, 'hvp': createRelationVaporizationEnthalpy}
    path = '/home/marina/PycharmProjects/combiff/wrkDir/tmp'

    def __init__(self, larsCode):
        self.larsCode = larsCode
        self.factories = {}

    def __str__(self):
        s = '\tLARSCODE: {}\n\tRELATIONS:'.format(self.larsCode)
        for rel in self.factories:
            s = '{} {},'.format(s, rel)
        s = s[:-1]
        return s

    def createFactory(larsCode, prop, rel, var, col):
        factory = dbsEntry.classes[prop](larsCode, prop, rel, var, col)
        return factory

    def addFactory(self, factory):
        rel = factory.rel
        if rel in self.factories:
            s = '{}: Relation already added.'.format(rel)
            print(s)
            sys.exit(666)
        self.factories[rel] = factory

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

                factory = dbsEntry.createFactory(larsCode, prop, rel, var, col)
                entry.addFactory(factory)

        return entry
    else:
        return obj