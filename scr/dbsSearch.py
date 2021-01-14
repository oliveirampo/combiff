from datetime import datetime
import configparser
import pandas as pd
import sys
import os


import myExceptions
import dbs


class dbsConfiguration():
    def __init__(self, fileName):
        if not os.path.exists(fileName):
            raise myExceptions.NoFile(fileName)

        config = configparser.ConfigParser()
        config.read(fileName)
        self.config = config

    def printSections(self):
        sections = self.config.sections()
        print(sections)

    def getPath(self):
        path = self.config.get('info', 'path')
        return path

    def getDbsFileName(self):
        fileName = self.config.get('info', 'dbsFileName')
        return fileName

    def getPropCod(self):
        # sections = self.config.sections()
        propCod = {}
        for option in self.config.options('relations'):
            propCod[option] = self.config.get('relations', option).split()
        return propCod

    def getDefaultPressure(self):
        val = self.config.get('defaultValues', 'pressure')
        val = float(val)
        return val

    def getDefaultTemperature(self):
        val = self.config.get('defaultValues', 'temperature')
        val = float(val)
        return val

    def getVariable(self, propName):
        try:
            var = self.config.get('globalVariables', propName)
        except configparser.NoOptionError:
            raise myExceptions.VariableNotDefined(propName)

        return var

    def getPropList(self):
        propList = self.config.get('plot', 'prop')
        propList = propList.split()
        return propList

    def getRelations(self, prop):
        relations = self.config.get('relations', prop)
        return relations

    def getNumberOfPoints(self):
        n = self.config.get('plot', 'nPoints')
        n = int(n)
        return n

    def isHvp(self, prop):
        propName = self.config.get('globalVariables', 'vaporizationEnthalpy')
        if prop == propName:
            return True
        else:
            return False

    def getOutFileName(self, typ):
        fileName = self.config.get('outputFiles', typ)
        return fileName


def run(dbsConfig):
    startTime = datetime.now()

    nArgs = len(sys.argv)
    if nArgs != 3:
        raise myExceptions.ArgError(nArgs, 3)

    # I decided to keep both files (molListFile and enuMolFile)
    # because I can easily comment out entries in molListFile
    # just be careful with "#" in smiles strings
    # 00_file.lst
    molListFile = sys.argv[2]
    # molList = np.genfromtxt(molListFile, dtype=None, encoding='utf-8')
    molList = pd.read_csv(molListFile, sep='\s+', header=None, names=['frm', 'cas', 'nam', 'inchi', 'smiles'])

    dbsFileName = dbsConfig.getDbsFileName()
    dbsEntries = dbs.getDbsEntries(dbsFileName)

    print('Getting Identifiers')
    path = dbsConfig.getPath()
    molList = getCids(dbsEntries, molList, path)
    # getCidsofSynonyms(dbsEntries, molList)

    print('Getting Data')
    propCod = dbsConfig.getPropCod()
    tables = getData(molList, dbsEntries, path, propCod)

    print('Writing Data')
    writeData(tables)

    print('\nTotal running time: {}\n'.format(datetime.now() - startTime))


def getCids(dbsEntries, molList, path):
    for larsCode in dbsEntries:
        # print(larsCode)
        factory = dbsEntries[larsCode].getFactory('cpd')
        parser = factory.createParser()
        dfTable = parser.readRelation(path)

        # TODO - match by smiles and inchikey (make it more general)
        # print('\tmatch by cas')
        molList = mergeByCas(larsCode, molList, dfTable)

        # print('\tmatch by name')
        mergeByName(larsCode, molList, dfTable)

    return molList


# function ignores if there are more than 1 match
def mergeByCas(larsCode, molList, dfTable):
    if not 'cas' in dfTable:
        return molList

    dfTable = dfTable[['cid', 'cas', 'nam']]
    dfTable.columns = ['cid_cas_' + larsCode, 'cas', 'nam_' + larsCode]
    dfTable = dfTable.loc[dfTable['cas'] != '%']

    df = pd.merge(molList, dfTable, how='left', on='cas')
    return df


def mergeByName(larsCode, molList, dfTable):
    dfTable = dfTable[['cid', 'nam']]
    dfTable.columns = ['cid_nam_' + larsCode, 'nam']
    dfTable = dfTable.loc[dfTable['nam'] != '%']

    df = pd.merge(molList, dfTable, how='left', on='nam')
    return df


def getCidsofSynonyms(dbsEntries, molList):
    factory = dbs.dbsEntry.getFactory('cpd')

    for larsCode in dbsEntries:
        compoundRelation = dbsEntries[larsCode].getCompoundRelation()
        parser = factory.createParser()
        dfTable = parser.readRelation(larsCode, compoundRelation)

        # match by name from cas match
        for idx, row in molList.iterrows():
            cid_cas = row['cid_cas_' + larsCode]
            cid_cas = str(cid_cas)

            if cid_cas == 'nan':
                cid_nam = row['cid_nam_' + larsCode]
                cid_nam = str(cid_nam)

                if cid_nam == 'nan':
                    print('TODO: Search for more identifiers using alternative names found.')
                    sys.exit(666)


def getData(molList, dbsEntries, path, propCod):
    columns = molList.columns.tolist()
    tables = {}

    for prop in propCod:
        tables[prop] = {}

        for larsCod in propCod[prop]:
            #print('\t', larsCod)
            tables[prop][larsCod] = []

            factory = dbsEntries[larsCod].getFactory(prop)
            parser = factory.createParser()
            dfTable = parser.readRelation(path)

            # match by cid
            # there will be lots of repeated data
            cid_options = [col for col in columns if col.startswith('cid_') and col.endswith(larsCod)]
            for cid in cid_options:
                tab = pd.merge(molList, dfTable, how='left', left_on=cid, right_on='cid')
                tab = tab.dropna(subset=[cid])
                tables[prop][larsCod].append(tab)

    return tables


def getFileName(prop, larsCode):
    fileName = 'data/{}_{}.csv'.format(prop, larsCode)
    return fileName


def writeData(tables):
    for prop in tables:
        #print(prop)
        for larsCode in tables[prop]:
            #print('\t', larsCode)

            tab = tables[prop][larsCode]
            if len(tab) == 0:
                continue

            df = pd.concat(tab, ignore_index=True)
            df = df.drop_duplicates(keep='last')

            df.reset_index(inplace=True, drop=True)

            fileName = getFileName(prop, larsCode)
            if df.shape[0] == 0:
                continue

            df.to_csv(fileName, index=False)


def readData(propCod):
    tables = {}

    for prop in propCod:
        tables[prop] = {}

        for larsCod in propCod[prop]:
            fileName = getFileName(prop, larsCod)

            try:
                tab = pd.read_csv(fileName)
            except FileNotFoundError:
                raise myExceptions.NoFile(fileName)

            tables[prop][larsCod] = tab

    return tables


