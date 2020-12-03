from collections import OrderedDict
from datetime import datetime
import pandas as pd
import numpy as np
import json
import sys
import os


import dbs
import myExceptions
from molecule import moleculeDecoder


#00_mol_A.lst ../../00_originals/mol_A.p


def run():
    startTime = datetime.now()

    nArgs = len(sys.argv)
    if nArgs != 4:
        raise myExceptions.ArgError(nArgs, 4)

    # I decided to keep both files (molListFile and enuMolFile)
    # because I can easily comment out entries in molListFile
    # 00_file.lst
    molListFile = sys.argv[2]
    # molList = np.genfromtxt(molListFile, dtype=None, encoding='utf-8')
    molList = pd.read_csv(molListFile, sep='\s+', header=None, names=['cod', 'frm', 'cas', 'nam', 'inchi', 'smiles'])

    # file.json
    enuMolFile = sys.argv[3]
    if not os.path.exists(enuMolFile):
        raise myExceptions.NoFile(enuMolFile)

    with open(enuMolFile) as jsonFile:
        enuData = json.load(jsonFile, object_hook=moleculeDecoder)

    data = OrderedDict()
    # check if 00_dbsData.json already exists
    # dbsDataFile = 'inp/00_dbsData.json'
    # if os.path.exists(dbsDataFile):
    #     with open(dbsDataFile) as jsonFile:
    #         data = json.load(jsonFile, object_hook=moleculeDecoder)

    # getData(molList, enuData, data, {})
    dbsFileName = 'inp/dbs.json'
    dbsEntries = getDbsEntries(dbsFileName)

    molList = getCids(dbsEntries, molList)
    # getCidsofSynonyms(dbsEntries, molList)

    propCod = {'dns': ['YA14.6'], 'hvp': ['AC16.1']}
    tables = getData(propCod, molList, dbsEntries)

    writeData(tables)

    print('\nTotal running time: {}\n'.format(datetime.now() - startTime))


def getDbsEntries(fileName):
    if not os.path.exists(fileName):
        raise myExceptions.NoFile(fileName)

    with open(fileName) as jsonFile:
        dbsEntries = json.load(jsonFile, object_hook=dbs.dbsEntryDecoder)
        # for entry in dbsEntries: print(dbsEntries[entry])

        return dbsEntries


def getCids(dbsEntries, molList):
    for larsCode in dbsEntries:
        factory = dbsEntries[larsCode].getFactory('cpd')
        parser = factory.createParser()
        dfTable = parser.readRelation()

        # match by cas
        df = dfTable[['cid', 'cas', 'nam']]
        df.columns = ['cid_cas_' + larsCode, 'cas', 'nam_' + larsCode]
        molList = pd.merge(molList, df, how='left', on='cas')

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


def getData(propCod, molList, dbsEntries):
    # print(molList)
    columns = molList.columns.tolist()
    tables = {}

    for prop in propCod:
        tables[prop] = {}

        for larsCod in propCod[prop]:
            tables[prop][larsCod] = []

            factory = dbsEntries[larsCod].getFactory(prop)
            parser = factory.createParser()
            dfTable = parser.readRelation()

            # match by cid
            # there will be lots of repeated data
            cid_options = [col for col in columns if col.startswith('cid_') and col.endswith(larsCod)]
            for cid in cid_options:
                tab = pd.merge(molList, dfTable, how='left', left_on=cid, right_on='cid')
                tab = tab.dropna(subset=[cid])
                tables[prop][larsCod].append(tab)

    return tables


def writeData(tables):
    dfAll = pd.DataFrame()

    for prop in tables:
        for larsCode in tables[prop]:

            tab = tables[prop][larsCode]
            df = pd.concat(tab, ignore_index=True)
            df = df.drop_duplicates(keep='last')

            fileName = 'data/{}_{}.csv'.format(prop, larsCode)
            df.to_csv(fileName)


