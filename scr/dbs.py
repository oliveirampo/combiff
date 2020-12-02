from abc import ABC, abstractmethod
import sys


import dbsIO


class abstractFactoryRelation(ABC):
    @abstractmethod
    def createRelation(self, rel, var, col):
        pass

    @abstractmethod
    def createParser(self):
        pass


class createRelationCompound(abstractFactoryRelation):
    def createRelation(self, rel, var, col):
        return dbsCompound(rel, var, col)

    def createParser(self):
        return dbsIO.parserDefault()


class createRelationDensity(abstractFactoryRelation):
    def createRelation(self, rel, var, col):
        return dbsDensity(rel, var, col)

    def createParser(self):
        return dbsIO.parserDefault()


class createRelationVaporizationEnthalpy(abstractFactoryRelation):
    def createRelation(self, rel, var, col):
        return dbsVaporizationEnthalpy(rel, var, col)

    def createParser(self):
        return dbsIO.parserDefault()


class dbsRelation(ABC):
    def __init__(self, rel, var, columns):
        self.rel = rel
        self.var = var
        self.columns = columns

    def __str__(self):
        s = '{} {}'.format(self.rel, self.var)
        return s


class dbsCompound(dbsRelation):
    def __init__(self, rel, var, columns):
        super(dbsCompound, self).__init__(rel, var, columns)


class dbsDensity(dbsRelation):
    def __init__(self, rel, var, columns):
        super(dbsDensity, self).__init__(rel, var, columns)


class dbsVaporizationEnthalpy(dbsRelation):
    def __init__(self, rel, var, columns):
        super(dbsVaporizationEnthalpy, self).__init__(rel, var, columns)


class dbsEntry():
    classes = {'cpd': createRelationCompound, 'dns': createRelationDensity, 'hvp': createRelationVaporizationEnthalpy}
    path = '/home/marina/PycharmProjects/combiff/wrkDir/tmp'

    def __init__(self, larsCode):
        self.larsCode = larsCode
        self.relations = {}

    def __str__(self):
        s = '\tLARSCODE: {}\n\tRELATIONS:'.format(self.larsCode)
        for rel in self.relations:
            s = '{} {},'.format(s, rel)
        s = s[:-1]
        return s

    def addRelation(self, relation):
        rel = relation.rel
        if rel in self.relations:
            s = '{}: Relation already added.'.format(rel)
            print(s)
            sys.exit(666)
        self.relations[rel] = relation

    def getCompoundRelation(self):
        return self.relations['cpd']


def dbsEntryDecoder(obj):
    if '__type__' in obj and obj['__type__'] == 'dbsEntry':
        entry = dbsEntry(obj['larsCode'])

        classes = dbsEntry.classes
        for prop in classes:
            if prop in obj:
                rel = obj[prop]['rel']
                var = obj[prop]['var']
                col = obj[prop]['col']

                factory = classes[prop]()
                relation = factory.createRelation(rel, var, col)
                entry.addRelation(relation)

        return entry
    else:
        return obj