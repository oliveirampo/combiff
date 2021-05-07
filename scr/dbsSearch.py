"""Get available data in DBS for a given molecule.

Classes:
    dbsConfiguration
Methods:
    run(dbsConfig)
    getCids(dbsEntries, molList, path)
    mergeByCas(larsCode, molList, dfTable)
    mergeByName(larsCode, molList, dfTable)
    getCidsofSynonyms(dbsEntries, molList)
    getData(molList, dbsEntries, path, propCod)
    getFileName(prop, larsCode)
    writeData(tables)
    readData(propCod)
"""


import configparser
import pandas as pd
import numpy as np
import sys
import os


import getIdentifiers
import myExceptions
import dbsSearch
import dbs


class dbsConfiguration:
    """DBS configuration object.

    Attributes:
        config: (configparser.ConfigParser) Configuration file parser.
    Methods:
        printSections()
        getPath()
        getDbsFileName()
        getPropCod()
        getDefaultPressure()
        getDefaultTemperature()
        getVariable(propName)
        getPropListToBePlotted()
        getPropListToBeWrittenToFile()
        getNumberOfPoints()
        isHvp(prop)
        getOutFileName(typ)
    """

    def __init__(self, fileName):
        """Constructs all the necessary attributes for the dbsConfiguration object.
        A configuration file consists of sections, lead by a "[section]" header,
        and followed by "name: value" entries, with continuations and such in
        the style of RFC 822.
        See 'configparser' documentation for more information.

        :param fileName: (str) File name from which configuration is read.
        """

        if not os.path.exists(fileName):
            raise myExceptions.NoFile(fileName)

        config = configparser.ConfigParser()
        config.read(fileName)
        self.config = config

    def printSections(self):
        """Prints a list of section names, excluding [DEFAULT]"""

        sections = self.config.sections()
        print(sections)

    def getPath(self):
        """Return path to tables with data from different sources.

        :return:
            path: (str) Path to tables with data from different sources.
        """

        path = self.config.get('info', 'path')
        return path

    def getDbsFileName(self):
        """Returns name of DNS json file (dbs.jon),
        which contains information about how to read tables in tmpDir/ directory.

        :return:
            fileName: (str) File name.
        """

        fileName = self.config.get('info', 'dbsFileName')
        return fileName

    def getPropCod(self):
        """Returns dictionary with
            key - propperty code
            value - associated LARS codes.

        :return:
            propCod: (dict) Dictionary that maps property code to LARDS codes.
        """

        # sections = self.config.sections()
        propCod = {}
        for option in self.config.options('relations'):
            propCod[option] = self.config.get('relations', option).split()
        return propCod

    def getDefaultPressure(self):
        """Returns default pressure defined in configuration file.

        :return:
            val: (float) Default pressure.
        """

        val = self.config.get('defaultValues', 'pressure')
        val = float(val)
        return val

    def getDefaultTemperature(self):
        """Returns default temperature defined in configuration file.

        :return:
            val: (str) Default temperature.
        """

        val = self.config.get('defaultValues', 'temperature')
        val = float(val)
        return val

    def getVariable(self, propName):
        """Returns code given property defined in globalVariables section.

        :param propName:
        :return:
            var: (str) Code.
        """

        try:
            var = self.config.get('globalVariables', propName)
        except configparser.NoOptionError:
            raise myExceptions.VariableNotDefined(propName)

        return var

    def getPropListToBePlotted(self):
        """Returns list of properties to be plotted.

        :return:
            propList: (list) List of properties.
        """

        propList = self.config.get('plot', 'prop')
        propList = propList.split()
        return propList

    def getPropListToBeWrittenToFile(self):
        """Returns list of properties to be plotted.

        :return:
            propList: (list) List of properties.
        """

        propList = self.config.get('writeProperty', 'prop')
        propList = propList.split()
        return propList

    def getNumberOfPoints(self):
        """Returns number of points to add in curve to be plotted.

        :return:
            n: (int) Number of points.
        """

        n = self.config.get('plot', 'nPoints')
        n = int(n)
        return n

    def isHvp(self, prop):
        """Returns if property if Vaporization Enthalpy.

        :param prop: (str) Property code.
        :return: (boolean)
        """

        propName = self.config.get('globalVariables', 'vaporizationEnthalpy')
        if prop == propName:
            return True
        else:
            return False

    def getOutFileName(self, typ):
        """Returns name of output file given type.

        :param typ: (str) Type of output file.
        :return:
            fileName: (str) Name out output file.
        """

        try:
            fileName = self.config.get('outputFiles', typ)
        except configparser.NoOptionError:
            raise myExceptions.VariableNotDefined(typ)
        return fileName

    def getFamilyCode(self):
        """Returns one letter family code."""

        try:
            fileName = self.config.get('writeProperty', 'family')
        except configparser.NoOptionError:
            raise myExceptions.VariableNotDefined('family')
        return fileName

    def getPlotDir(self):
        """Returns directory where plots will be saved."""

        plotDir = self.config.get('directories', 'plotDir')
        return plotDir

    def getSimDir(self):
        """Returns directory with simulation results."""

        simDir = self.config.get('directories', 'simulationDir')
        return simDir


def read_mol_list(file_name):
    molList = pd.read_csv(file_name, sep='\s+', header=None, names=['frm', 'cas', 'nam', 'inchi', 'smiles'])
    return molList


def run(dbsConfig):
    """Get identifiers and all available data for a given molecule in DBS
    and save data to data/ directory.

    :param dbsConfig: (dbsConfiguration object) DBS configuration object.
    """

    nArgs = len(sys.argv)
    if nArgs != 3:
        raise myExceptions.ArgError(nArgs, 3)

    # I decided to keep both files (molListFile and enuMolFile)
    # because I can easily comment out entries in molListFile,
    # because of "#" in smiles strings I cannot use this symbol to start a comment.

    # 00_file.lst
    molListFile = sys.argv[2]
    # molList = np.genfromtxt(molListFile, dtype=None, encoding='utf-8')
    if not os.path.exists(molListFile):
        raise myExceptions.NoFile(molListFile)
    molList = read_mol_list(molListFile)

    dbsFileName = dbsConfig.getDbsFileName()
    dbsEntries = dbs.getDbsEntries(dbsFileName)

    print('Getting Identifiers')
    path = dbsConfig.getPath()
    molList = getCids(dbsEntries, molList, path)

    getCidsofSynonyms(dbsConfig, dbsEntries, molList, path)

    print('Getting Data')
    propCod = dbsConfig.getPropCod()
    tables = getData(dbsEntries, molList, path, propCod)

    print('Writing Data')
    writeData(tables)


def getCids(dbsEntries, molList, path):
    """Return augmented list of molecules with matched identifiers (CIDs) from DBS.

    :param dbsEntries: (dict) DBS entries. One for each LARS code.
    :param molList: (pandas DataFrame) Table with molecules and general identifiers (formula, cas, name, inchi, smiles).
    :param path: (str) Path to directory with DBS tables (tmp/).
    :return:
        molList: (pandas DataFrame) Table with molecules all identifiers (general ones + DBS ones).
    """

    for larsCode in dbsEntries:
        # print(larsCode)
        factory = dbsEntries[larsCode].getFactory('cpd')
        parser = factory.getParser()
        dfTable = parser.readRelation(path)

        # TODO - match by smiles and inchikey (make it more general)
        # print('\tmatch by cas')
        molList = mergeByCas(larsCode, molList, dfTable)

        # print('\tmatch by name')
        mergeByName(larsCode, molList, dfTable)

    return molList


def mergeByCas(larsCode, molList, dfTable):
    """Merges two tables by the CAS number.
    1st table - List of molecules and identifiers extracted from PubChem.
    2ns table - List of molecules and CID from DBS for the associated LARS code.

    :param larsCode: (str) LARS code.
    :param molList: (pandas DataFrame) Table with molecules and identifiers.
    :param dfTable: (pandas DataFrame) Compound table associated with the given LARS code.
    :return:
        df: (pandas DataFrame) Augmented table with molecules and identifiers.

    Function ignores rows if there is more than 1 match.
    """

    if 'cas' not in dfTable:
        return molList

    dfTable = dfTable[['cid', 'cas', 'nam']]
    dfTable.columns = ['cid_cas_' + larsCode, 'cas', 'nam_' + larsCode]
    dfTable = dfTable.loc[dfTable['cas'] != '%']

    df = pd.merge(molList, dfTable, how='left', on='cas')

    if df.shape[0] != molList.shape[0]:
        print('Multiple rows with same cas:', larsCode)
        res = df[df.isin(df[df.duplicated(subset=['cas'])])].dropna(subset=['cas'])
        print(np.unique(res['cas'].values))
        sys.exit(1)

    return df


def mergeByName(larsCode, molList, dfTable):
    """Merges two tables by the name.
    1st table - List of molecules and identifiers extracted from PubChem.
    2ns table - List of molecules and CID from DBS for the associated LARS code.

    :param larsCode: (str) LARS code.
    :param molList: (pandas DataFrame) Table with molecules and identifiers.
    :param dfTable: (pandas DataFrame) Compound table associated with the given LARS code.
    :return:
        df: (pandas DataFrame) Augmented table with molecules and identifiers.

    Function ignores rows if there is more than 1 match.
    """

    dfTable = dfTable[['cid', 'nam']]
    dfTable.columns = ['cid_nam_' + larsCode, 'nam']
    dfTable = dfTable.loc[dfTable['nam'] != '%']

    df = pd.merge(molList, dfTable, how='left', on='nam')
    return df


def getCidsofSynonyms(dbsConfig, dbsEntries, molList, path):
    """Finds synonyms and adds their CIDs from DBS.

    :param dbsEntries: (dict) DBS entries. One for each LARS code.
    :param molList: (pandas DataFrame) Table with molecules and identifiers.
    """

    propCod = dbsConfig.getPropCod()
    relations = get_relations(propCod)

    for rel in molList.columns:
        if rel.startswith('nam_'):
            larsCode = rel.split('_')[-1]

            factory = dbsEntries[larsCode].getFactory('cpd')
            parser = factory.getParser()
            dfTable = parser.readRelation(path)

            dfTable = dfTable[dfTable['nam'] != '%'][['cid', 'nam']]
            dfTable.columns = ['cid', 'synonym']

            molList = pd.merge(molList, dfTable, left_on='nam', right_on='synonym', how='left')

            result = molList[molList['cid'].notnull()]

            if result.shape[0] == 0:
                molList = molList.drop(columns=['cid', 'synonym'], axis=1)

            else:
                for idx, row in molList[molList['cid'].notnull()].iterrows():
                    cid = row['cid']
                    nam_col = 'cid_nam_' + larsCode
                    cas_col = 'cid_cas_' + larsCode

                    # Check if cid is already in table
                    # if not add to new column
                    saved_cid = molList.loc[idx, cas_col]
                    if saved_cid != cid:
                        if nam_col not in molList.columns:
                            molList[nam_col] = ''
                            molList.loc[idx, nam_col] = cid
                            print('\t', cid)
                        else:
                            sys.exit('TODO: add matches.')

                molList = molList.drop(columns=['cid', 'synonym'], axis=1)

    # sys.exit('STOP')


def getData(dbsEntries, molList, path, propCod):
    """

    :param dbsEntries: (dict) DBS entries. One for each LARS code.
    :param molList: (pandas DataFrame) Table with molecules and identifiers.
    :param path: (str) Path to directory with DBS tables (tmp/).
    :param propCod: (str) Property code.
    :return:
        tables: (dict) Tables of property and source (LARS code).
                {property code {LARS code: [values]}}
    """

    columns = molList.columns.tolist()
    tables = {}

    for prop in propCod:
        tables[prop] = {}

        for larsCod in propCod[prop]:
            # print('\t', larsCod)
            tables[prop][larsCod] = []

            factory = dbsEntries[larsCod].getFactory(prop)
            parser = factory.getParser()
            dfTable = parser.readRelation(path)

            # match by cid
            # there will be lots of repeated data
            cid_options = [col for col in columns if col.startswith('cid_') and col.endswith(larsCod)]
            for cid in cid_options:

                # remove unnecessary columns
                keep_columns = ['frm', 'cas', 'nam', 'inchi', 'smiles']
                for col in columns:
                    if larsCod in col:
                        keep_columns.append(col)

                tab = pd.merge(molList[keep_columns], dfTable, how='left', left_on=cid, right_on='cid')
                tab = tab.dropna(subset=[cid])

                tab = parser.removeEmptyValues(factory.var, tab)

                tables[prop][larsCod].append(tab)

    return tables


def getFileName(prop, larsCode):
    """Returns name of file given property and LARS code where data will be saved to.

    :param prop: (str) Property code.
    :param larsCode: (str) LARS code.
    :return:
        fileName: (str) Name of file.
    """

    fileName = 'data/{}_{}.csv'.format(prop, larsCode)
    return fileName


def writeData(tables):
    """Writes data extracted from DBS tables to files.

    :param tables: (dict) Tables of property and source (LARS code).
                    {property code {LARS code: [values]}}
    """

    for prop in tables:
        # print(prop)
        for larsCode in tables[prop]:
            # print('\t', larsCode)

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
    """Reads files with data extracted from DBS.

    :param propCod: (str) Property code.
    :return:
        tables: (dict) Tables of property and source (LARS code).
                    {property code {LARS code: [values]}}
    """

    tables = {}

    for prop in propCod:
        tables[prop] = {}

        for larsCod in propCod[prop]:
            fileName = getFileName(prop, larsCod)

            if not os.path.exists(fileName):
                continue

            try:
                tab = pd.read_csv(fileName)
            except FileNotFoundError:
                raise myExceptions.NoFile(fileName)

            tables[prop][larsCod] = tab

    return tables


def runAlternative(dbsConfig, alternative_directory):
    molListFile = sys.argv[2]
    if not os.path.exists(molListFile):
        raise myExceptions.NoFile(molListFile)
    molList = read_mol_list(molListFile)

    dbsFileName = dbsConfig.getDbsFileName()
    dbsEntries = dbs.getDbsEntries(dbsFileName)

    print('Getting Identifiers')
    path = dbsConfig.getPath()
    molList = getCids(dbsEntries, molList, path)

    propCod = dbsConfig.getPropCod()
    relations = get_relations(propCod)

    for larsCode in relations:
        file_name = '../{}/{}_cpd.rel'.format(alternative_directory, larsCode)
        if os.path.exists(file_name):
            # print(file_name)
            alternative_data = pd.read_csv(file_name, sep='\s+')

            # Merge by smiles
            cid_smiles_col = 'cid_smiles_{}'.format(larsCode)
            current_table = alternative_data.loc[alternative_data['smiles'] != '%'][['cid', 'smiles']]
            current_table.columns = [cid_smiles_col, 'smiles']
            molList = pd.merge(molList, current_table, how='left', on='smiles')

            molList[cid_smiles_col] = molList[cid_smiles_col].str.lower()

            # check for duplicated columns
            if cid_smiles_col + '_x' in molList.columns:
                sys.exit('TODO: Duplicated columns for smiles')

            # Merge by CAS
            cid_cas_col = 'cid_cas_{}'.format(larsCode)
            current_table = alternative_data.loc[alternative_data['cas'] != '%'][['cid', 'cas']]
            current_table.columns = [cid_cas_col, 'cas']
            molList = pd.merge(molList, current_table, how='left', on='cas')

            # check for duplicated columns
            if cid_cas_col + '_x' in molList.columns:
                df_tmp = molList[molList[cid_cas_col + '_y'].notnull()]

                if df_tmp.shape[0] == 0:
                    molList = molList.drop(columns=[cid_cas_col + '_y'], axis=1)

                    columns = molList.columns.values
                    for i, col in enumerate(columns):
                        if col == cid_cas_col + '_x':
                            columns[i] = cid_cas_col

                    molList.columns = columns

                else:
                    # print(file_name)
                    # print(molList[[cid_cas_col + '_x', cid_cas_col + '_y']].head(2))
                    # print(df_tmp[[cid_cas_col + '_x', cid_cas_col + '_y']].head(2))
                    sys.exit('TODO: Duplicated columns for CAS')

            # sys.exit('STOP')

    print('Getting Data')
    propCod = dbsConfig.getPropCod()
    tables = getData(dbsEntries, molList, path, propCod)

    print('Writing Data')
    writeData(tables)


def get_relations(propCod):
    relations = []
    for prop in propCod:
        for rel in propCod[prop]:
            relations.append(rel)
    relations = list(set(relations))
    return relations
