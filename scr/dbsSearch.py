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

    def getProps(self):
        props = self.config.get('plot', 'prop')
        props = props.split()
        return props

    def getVariable(self, propName):
        var = self.config.get('globalVariables', propName)
        return var

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
    molList = pd.read_csv(molListFile, sep='\s+', header=None, names=['cod', 'frm', 'cas', 'nam', 'inchi', 'smiles'])

    dbsFileName = dbsConfig.getDbsFileName()
    dbsEntries = dbs.getDbsEntries(dbsFileName)

    path = dbsConfig.getPath()
    molList = getCids(dbsEntries, molList, path)
    # getCidsofSynonyms(dbsEntries, molList)

    propCod = dbsConfig.getPropCod()
    tables = getData(molList, dbsEntries, path, propCod)

    writeData(tables)

    print('\nTotal running time: {}\n'.format(datetime.now() - startTime))


def getCids(dbsEntries, molList, path):
    for larsCode in dbsEntries:
        factory = dbsEntries[larsCode].getFactory('cpd')
        parser = factory.createParser()
        dfTable = parser.readRelation(path)

        # TODO - match my smiles and inchikey (make it more general)
        try:
            # match by cas
            df = dfTable[['cid', 'cas', 'nam']]
            df.columns = ['cid_cas_' + larsCode, 'cas', 'nam_' + larsCode]
            molList = pd.merge(molList, df, how='left', on='cas')
        except KeyError as err:
            # print(larsCode, err)
            pass

        # match by name
        df = dfTable[['cid', 'nam']]
        df.columns = ['cid_nam_' + larsCode, 'nam']
        molList = pd.merge(molList, df, how='left', on='nam')

    return molList


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
    # print(molList)
    columns = molList.columns.tolist()
    tables = {}

    # propCod = dbs.dbsEntry.propCod
    for prop in propCod:
        tables[prop] = {}

        for larsCod in propCod[prop]:
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
        for larsCode in tables[prop]:

            tab = tables[prop][larsCode]
            df = pd.concat(tab, ignore_index=True)
            df = df.drop_duplicates(keep='last')

            df.reset_index(inplace=True, drop=True)

            fileName = getFileName(prop, larsCode)
            df.to_csv(fileName)


def readData(propCod):
    tables = {}

    for prop in propCod:
        tables[prop] = {}

        for larsCod in propCod[prop]:
            fileName = getFileName(prop, larsCod)
            tab = pd.read_csv(fileName, index_col=0)

            tables[prop][larsCod] = tab

    return tables


