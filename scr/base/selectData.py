"""Module used to plot all data available and select data points manually.

Methods:
    manualSelection(dbsConfig)
    getDbsData(prop, larsCode, dbsEntries)
    plotSelectedData(propVar, selectedData, defaultPressure)
    getTransitionPoint(propName, defaultPressure, dbsConfig, dbsEntries, propCod, tables, smiles)
    addVaporPressure(propCod, blp, tables, smiles, dbsEntries, x_values, pre_values, src_values)
    getFirstValueOfTransitionPoint(data)
    getHighestTransitionPoint(var, data)
    convertArrays(prop, x_values, y_values, pre_values, src_values, fid_values, met_values)
    getSelectedDataForProperty(dbsConfig)
    getSelectedData(dbsConfig)
    addProperty(prop, dfTmp, selectedData)
    assignCode(dbsConfig)
"""

from rdkit.Chem import rdMolDescriptors
import matplotlib.pyplot as plt
from rdkit import Chem
import pandas as pd
import numpy as np
import sys
import os

from scr.base.dbsRelation import getUnit
from scr.base import myDataStructure
from scr.base import myExceptions
from scr.base import family_utils
from scr.base import dbsSearch
from scr.base import plotData
from scr.base import utils
from scr.base import dbs
from scr.base import IO


def manualSelection(dbsConfig):
    """
    Plots properties for each molecule as function of temperature,
    and writes manually selected data points to file.

    :param dbsConfig: (dbsConfiguration object) DBS configuration object.
    """

    nArgs = len(sys.argv)
    if nArgs != 3:
        raise myExceptions.ArgError(nArgs, 3)

    # 00_file.lst
    molListFile = sys.argv[2]
    molList = pd.read_csv(molListFile, sep='\s+', header=None, names=['frm', 'cas', 'nam', 'inchi', 'smiles'])

    dbsFileName = dbsConfig.getDbsFileName()
    dbsEntries = dbs.getDbsEntries(dbsFileName)

    # get tables with data
    properties = dbsConfig.getPropListToBePlotted()
    propCod = dbsConfig.getPropCod()
    tables = dbsSearch.readData(propCod)

    defaultPressure = dbsConfig.getDefaultPressure()
    tem_room = dbsConfig.getDefaultTemperature()
    nPoints = dbsConfig.getNumberOfPoints()

    allSelectedData = IO.openSelectedDataFile(dbsConfig)

    plotValues(dbsConfig, dbsEntries, defaultPressure, tem_room, nPoints, properties, propCod, tables, molList,
             allSelectedData, True, True, False, False)


def plotValues(dbsConfig, dbsEntries, defaultPressure, tem_room, nPoints, properties, propCod, tables, molList,
             allSelectedData, isSelectData, isWriteData, isSavePlot, isPlotSim):
    """Plots all data.

    :param dbsConfig: (dbsConfiguration object) DBS configuration object.
    :param dbsEntries: (dict of dbsEntry objects) Dictionary that maps LARS code to dbsEntry objects.
    :param defaultPressure: (float) Default pressure.
    :param tem_room: (float) Room temperature.
    :param nPoints: (int) Number of points to plot.
    :param properties: (dict) Dictionary of properties.
    :param propCod: (dict) Map from property code to LARS code available for this property.
    :param tables: tables: (dict) Dictionary of pandas DataFrame that maps property codes to data.
    :param molList: (pandas DataFrame) Table with molecules and identifiers.
    :param allSelectedData: (dict of selectedData) Dict that maps SMILES to object that contains the selected data.
    :param isSelectData: (float) Allow data selection.
    """

    for prop in properties:
        propUnit = getUnit(prop)

        smilesVariable = dbsConfig.getVariable('smiles')

        simDir = dbsConfig.getSimDir()
        simDataFile = '{}/{}.dat'.format(simDir, prop)
        if os.path.exists(simDataFile):
            allSimData = pd.read_csv(simDataFile, sep='\s+', names=[smilesVariable, 'tem', 'val'])

        for idx, row in molList.iterrows():
            smiles = row[smilesVariable]

            selectedData = myDataStructure.SelectedData(smiles)
            if smiles in allSelectedData:
                selectedData = allSelectedData[smiles]

            mlpVar, mlp = getTransitionPoint('meltingPoint', defaultPressure, dbsConfig, dbsEntries, propCod, tables,
                                             smiles)
            blpVar, blp = getTransitionPoint('boilingPoint', defaultPressure, dbsConfig, dbsEntries, propCod, tables,
                                             smiles)
            tem_cri_var, tem_cri = getTransitionPoint('criticalTemperature', defaultPressure, dbsConfig, dbsEntries,
                                                      propCod, tables, smiles)

            selectedData.addMeltingPoint(mlp)
            selectedData.addBoilingPoint(blp)
            selectedData.addCriticalTemperature(tem_cri)

            fig, ax = plotData.plotStart(mlpVar, mlp, blpVar, blp, tem_cri_var, tem_cri)
            if selectedData.hasProperty(prop):
                plotSelectedData(prop, selectedData, defaultPressure)

            x_values = []
            y_values = []
            pre_values = []
            src_values = []
            fid_values = []
            met_values = []

            if dbsConfig.isHvp(prop):
                plotData.plotHvb(propCod, tables, smiles, dbsEntries, defaultPressure, x_values, y_values, src_values,
                                 pre_values, fid_values, met_values)

            for larsCode in propCod[prop]:
                if larsCode not in tables[prop]:
                    continue

                tab = tables[prop][larsCode]
                tab = tab.loc[tab['smiles'] == smiles]

                if tab.shape[0] == 0:
                    continue

                dbsRelation, dbsFactory = getDbsData(prop, larsCode, dbsEntries)
                marker = dbsRelation.getMarker()
                color = dbsRelation.getColor()
                if prop in tab:
                    data = dbsRelation.getData(tab, defaultPressure)

                    pre, tem, val, fid, met = dbsRelation.getValues(data)

                    plotData.plotPoint(tem, val, larsCode, marker, color, linewidth=0.0)

                    for i in range(tem.shape[0]):
                        t = tem[i]
                        v = val[i]
                        f = fid[i]
                        m = met[i]

                        if f == '%':
                            f = ''
                        if m == '%':
                            m = ''

                        plotData.addCode(t, v, f)
                        plotData.addCode(t, v, m)

                    plotData.saveDataPts(tem, val, larsCode, fid, met, x_values, y_values, src_values, fid_values, met_values)
                    for p in pre: pre_values.append(p)

                equations = dbsFactory.createEquations(tab)
                for eq in equations:
                    X, Y, fid = eq.getData(tem_room, nPoints, dbsRelation.tem_convert)

                    if len(X) > 0:
                        fid = X.shape[0] * [fid[0]]
                        met = X.shape[0] * ['']

                        plotData.plotPoint(X, Y, larsCode, marker, color, linewidth=1.0)
                        plotData.saveDataPts(X, Y, larsCode, fid, met, x_values, y_values, src_values, fid_values,
                                             met_values)

                        f = fid[0]
                        if f == '%':
                            f = ''
                        plotData.addCode(X[0], Y[0], f)

                        for x in X: pre_values.append(defaultPressure)

            x_values, y_values, pre_values, src_values, fid_values, met_values = \
                convertArrays(prop, x_values, y_values, pre_values, src_values, fid_values, met_values)

            plotData.plotDetails(fig, row, propUnit)

            addVaporPressure(propCod, blp, tables, smiles, dbsEntries, x_values, pre_values, src_values)

            plotData.markPressure(x_values, y_values, pre_values, defaultPressure)

            # Plot simulated data
            if isPlotSim:
                print(smiles)
                simData = allSimData.loc[allSimData[smilesVariable] == smiles]
                if simData.shape[0] != 0:
                    for idxSim, rowSim in simData.iterrows():
                        simTem = rowSim['tem']
                        simVal = rowSim['val']
                        plotData.plotPoint(simTem, simVal, 'Simulated', 'd', 'red', linewidth=0.0)

            selectedPointsMask = []
            if len(x_values) != 0:
                selectedPointsMask = plotData.select_data(fig, ax, x_values, y_values)

            if isSelectData:
                selectedData.appendProperty(prop, selectedPointsMask, x_values, y_values, pre_values, src_values,
                                            fid_values, met_values)

            if isWriteData:
                if selectedData.hasData(prop):
                    print(smiles)
                    allSelectedData[smiles] = selectedData
                    IO.writeSelectedDataToJson(dbsConfig, allSelectedData)

            if isSavePlot:
                plotDir = dbsConfig.getPlotDir()
                code = selectedData.code

                if not code:
                    print('\nMolecule without code.\nRun option -assignCode for all molecules.')
                    sys.exit(1)

                plotFileName = '{}/plot_{}_{}.png'.format(plotDir, prop, code)
                fig.savefig(plotFileName)

            plt.close(fig)


def getDbsData(prop, larsCode, dbsEntries):
    """Returns dbs relation objects and factory object that creates equation objects.

    :param prop: (str) Property code.
    :param larsCode: (str) LARS code.
    :param dbsEntries: (dict of dbsEntry objects) Dictionary that maps LARS code to dbsEntry objects.
    :return:
        dbsRelation: (dbsRelation object)
        dbsFactory: (abstractFactoryRelation object) Object that creates relations and equations.
    """

    dbsEntry = dbsEntries[larsCode]
    dbsFactory = dbsEntry.getFactory(prop)
    dbsRelation = dbsFactory.getRelation()

    return dbsRelation, dbsFactory


def plotSelectedData(propVar, selectedData, defaultPressure):
    """Plots selected data.

    :param propVar: (str) Property code.
    :param selectedData: (selectedData object) Object that contains data selected.
    :param defaultPressure: (float) Vapor pressure.
    :return:
    """

    prop = selectedData.getProperty(propVar)
    pre = np.asarray(prop.pre, dtype=np.float32)
    tem = np.asarray(prop.tem, dtype=np.float32)
    val = np.asarray(prop.val, dtype=np.float32)
    src = np.asarray(prop.src)

    for i in range(val.shape[0]):
        plotData.plotPoint(tem[i], val[i], src[i], 'd', 'gold', linewidth=0.0)
    plotData.markPressure(tem, val, pre, defaultPressure)


def getTransitionPoint(propName, defaultPressure, dbsConfig, dbsEntries, propCod, tables, smiles):
    """Returns transition point.

    :param propName: (str) Property code: meltingPoint, boilingPoint or criticalTemperature.
    :param defaultPressure: (float) Default pressure.
    :param dbsConfig: (dbsConfiguration object) DBS configuration object.
    :param dbsEntries: (dict of dbsEntry objects) Dictionary that maps LARS code to dbsEntry objects.
    :param propCod: (dict) Map from property code to LARS code available for this property.
    :param tables: (dict) Dictionary of pandas DataFrame that maps property codes to data.
    :param smiles: (str) SMILES string.
    :return:
        prop: (str) Property code.
        data: (pandas DataFrame) Table with two columns - temperature and larsCode.
    """

    prop = dbsConfig.getVariable(propName)

    data = []
    for larsCode in propCod[prop]:
        if larsCode in tables[prop]:
            tab = tables[prop][larsCode]
            tab = tab.loc[tab['smiles'] == smiles]

            if tab.shape[0] != 0:
                dbsRelation, dbsFactory = getDbsData(prop, larsCode, dbsEntries)
                dat = dbsRelation.getData(tab, defaultPressure)

                values = dat.shape[0] * [[0.0, larsCode]]

                for i in range(dat.shape[0]):
                    values[i][0] = dat[i, dbsRelation.propCol]
                data.extend(values)

    data = pd.DataFrame(data, columns=['tem', 'larsCode'])
    return prop, data


def addVaporPressure(propCod, blp, tables, smiles, dbsEntries, x_values, pre_values, src_values):
    """

    :param propCod: (dict) Map from property code to LARS code available for this property.
    :param blp: (pandas DataFrame) Table with two columns [boiling point | LARS code].
    :param tables: (dict) Dictionary of pandas DataFrame that maps property codes to data.
    :param smiles: (str) SMILES string.
    :param dbsEntries: (dict of dbsEntry objects) Dictionary that maps LARS code to dbsEntry objects.
    :param x_values: (numpy.ndarray) Array with x values.
    :param pre_values: (numpy.ndarray) Vapor pressure values.
    :param src_values: (numpy.ndarray) Array with LARS code for each pressure value.
    :return:
    """

    prop = 'pvp'
    if prop not in propCod:
        return

    maxBlp = getHighestTransitionPoint('tem', blp)

    for larsCode in propCod[prop]:
        if larsCode not in tables[prop]:
            continue

        tab = tables[prop][larsCode]
        tab = tab.loc[tab['smiles'] == smiles]
        if tab.shape[0] == 0:
            continue

        dbsEntry = dbsEntries[larsCode]
        dbsFactory = dbsEntry.getFactory(prop)

        equations = dbsFactory.createEquations(tab)

        for eq in equations:
            for i in range(len(src_values)):
                src = src_values[i]

                if src == 'YA14.6':
                    tem = x_values[i]

                    if tem > maxBlp:
                        val = eq.getDataAt(tem)
                        if len(val) == 1:
                            val = val[0]
                            pre_values[i] = val


def getFirstValueOfTransitionPoint(data):
    """Returns first value of a given transition point.

    :param data: (pandas DataFrame) Table with temperature of transition point and LARS code.
    :return:
        (array) Array with first temperature and LARS code pair available.
    """

    if data.shape[0] != 0:
        row = data.iloc[0]
        tem = row['tem']
        src = row['larsCode']
        return [tem, src]
    return []


def getHighestTransitionPoint(var, data):
    """Returns maximum value in a list for a given transition point.

    :param var: (str) Column name in 'data' argument.
    :param data: (pandas DataFrame) Table with temperature of transition point and LARS code.
    :return:
    """

    maxVal = data[var].max()
    maxVal = float(maxVal)
    return maxVal


def convertArrays(prop, x_values, y_values, pre_values, src_values, fid_values, met_values):
    """Converts lists to numpy arrays.

    :param prop: (str) Property code.
    :param x_values: (list) X values.
    :param y_values: (list) Y values.
    :param pre_values: (list) Pressure value.
    :param src_values: (list) List of LARS codes.
    :param fid_values: (list) List of annotations if provided by the source.
    :param met_values: (list) List of methods if provided by the source.
    :return:
        x_values: (numpy array)
        y_values: (numpy array)
        pre_values: (numpy array)
        src_values: (numpy array)
        fid_values: (numpy array)
        met_values: (numpy array)
    """

    x_values = np.asarray(x_values)
    y_values = np.asarray(y_values)
    pre_values = np.asarray(pre_values)
    src_values = np.asarray(src_values)
    fid_values = np.asarray(fid_values)
    met_values = np.asarray(met_values)

    if len(x_values) != len(y_values):
        print('len(tem) != len({})\n'.format(prop))
        sys.exit(1)
    if len(x_values) != len(pre_values):
        print('len(tem) != len(pre) for prop = {}\n'.format(prop))
        sys.exit(1)

    return x_values, y_values, pre_values, src_values, fid_values, met_values


def getSelectedDataForProperty(dbsConfig):
    """Gets selected data from old format file
    and adds to already available output file.

    :param dbsConfig: (dbsConfiguration object) DBS configuration object.
    """

    fieFile = sys.argv[2]
    isomers = IO.readFieFile(fieFile)

    oldMolListFileName = sys.argv[3]    # 01_mol.lst

    molList = pd.read_csv(oldMolListFileName, sep='\s+',
              names=['cod', 'frm', 'cas', 'nam', 'inchi', 'smiles', 'old_cod', 'x', 'y', 'z'],
              usecols=['cod', 'inchi', 'smiles'])

    propCod = dbsConfig.getPropCod()
    allSelectedData = IO.openSelectedDataFile(dbsConfig)

    for prop in propCod:
        if prop == 'diffus':
            oldSelectedDataFileName = 'old_files/{}.src'.format('self_dif')
        elif prop == 'etd':
            oldSelectedDataFileName = 'old_files/{}.src'.format('vis_d')
        else:
            oldSelectedDataFileName = 'old_files/{}.src'.format(prop)

        if not os.path.exists(oldSelectedDataFileName):
            continue

        print(prop)

        srcData = pd.read_csv(oldSelectedDataFileName, sep='\s+', names=['cod', 'pre', 'tem', 'val', 'src'])
        srcData = srcData.loc[srcData['val'].astype(str) != '%']

        mergedData = pd.merge(srcData, molList, how='inner', on='cod')

        # get smiles in common
        isomers = utils.canonicalizeSmiles(isomers)
        mergedData = utils.canonicalizeSmiles(mergedData)
        mergedData = pd.merge(isomers, mergedData, how='inner', on='smiles')

        for smiles in allSelectedData:
            selectedData = allSelectedData[smiles]

            df = mergedData.loc[mergedData['smiles'] == smiles]

            nRows = df.shape[0]
            if nRows == 0:
                continue
            elif nRows != 1:
                print(df)
                print('More than one row for {} {}'.format(smiles, prop))
                sys.exit(1)

            pre = df['pre'].values
            tem = df['tem'].values
            val = df['val'].values
            src = df['src'].values

            fid = val.shape[0] * ['']
            met = val.shape[0] * ['']

            if not selectedData.hasProperty(prop):
                property = myDataStructure.Property(prop, pre, tem, val, src, fid, met)
                selectedData.addProperty(prop, property)

    IO.writeSelectedDataToJson(dbsConfig, allSelectedData)


def getSelectedData(dbsConfig):
    """Gets selected data from old format file for dns and hvp properties.

    :param dbsConfig:  (dbsConfiguration object) DBS configuration object.
    :return:
    """

    fieFile = sys.argv[2]
    isomers = IO.readFieFile(fieFile)

    oldMolListFileName = sys.argv[3]    # 01_mol.lst
    oldSelectedDataFileName = 'old_files/mol.src'

    molList = pd.read_csv(oldMolListFileName, sep='\s+',
                          names=['cod', 'frm', 'cas', 'nam', 'inchi', 'smiles', 'old_cod', 'x', 'y', 'z'],
                          usecols=['cod', 'inchi', 'smiles'])

    srcData = pd.read_csv(oldSelectedDataFileName, sep='\s+',
                          names=['cod', 'nam', 'cas', 'pre', 'tem', 'dns', 'dns_src', 'hvp', 'hvp_src',
                                 'mlp', 'mlp_src', 'blp', 'blp_src', 'tem_cri', 'tem_cri_src',
                                 'eps', 'eps_src', 'x', 'y', 'z'])

    srcData['letter'] = srcData['cod'].str[-1]
    srcData['cod'] = srcData['cod'].str[:-1]

    # remove row which start with \#
    for idx, row in srcData.iterrows():
        cod = row['cod']
        if cod[0] == '#':
            srcData.drop([idx], inplace=True)

    df = pd.merge(srcData, molList, how='left', on='cod')

    # get smiles in common
    isomers = utils.canonicalizeSmiles(isomers)
    df = utils.canonicalizeSmiles(df)
    df = pd.merge(isomers, df, how='left', on='smiles')

    # create selectedData
    # write to json file
    allSelectedData = {}
    codes = df['cod'].unique()
    for cod in codes:
        appendData = False
        dfTmp = df.loc[df['cod'] == cod]

        smiles = dfTmp['smiles'].unique()
        if smiles.shape[0] != 1:
            continue
        smiles = smiles[0]
        # print('COD:', cod)

        selectedData = myDataStructure.SelectedData(smiles)

        addProperty('mlp', dfTmp, selectedData)
        addProperty('blp', dfTmp, selectedData)
        addProperty('tem_cri', dfTmp, selectedData)
        addProperty('eps', dfTmp, selectedData)

        data = dfTmp[dfTmp['dns'].astype(str) != '%']
        data = data[data['dns'].astype(str) != '-']
        if data.shape[0] != 0:
            pre = data['pre'].values * 100.00
            tem = data['tem'].values
            dns = data['dns'].values
            src = data['dns_src'].values
            appendData = True

            fid = dns.shape[0] * ['']
            met = dns.shape[0] * ['']

            property = myDataStructure.Property('dns', pre, tem, dns, src, fid, met)
            selectedData.addProperty('dns', property)

        data = dfTmp[dfTmp['hvp'].astype(str) != '%']
        data = data[data['hvp'].astype(str) != '-']
        if data.shape[0] != 0:
            pre = data['pre'].values * 100.00
            tem = data['tem'].values
            hvp = data['hvp'].values
            src = data['hvp_src'].values
            appendData = True

            fid = hvp.shape[0] * ['']
            met = hvp.shape[0] * ['']

            property = myDataStructure.Property('hvp', pre, tem, hvp, src, fid, met)
            selectedData.addProperty('hvp', property)

        if appendData:
            allSelectedData[smiles] = selectedData
            IO.writeSelectedDataToJson(dbsConfig, allSelectedData)


def addProperty(prop, dfTmp, selectedData):
    """Adds property to selected data object.

    :param prop: (str) Property code.
    :param dfTmp: (pandas DataFrame) Table with data.
    :param selectedData: selectedData: (selectedData object) Object that contains data selected.
    :return:
    """

    val = dfTmp[prop].unique()
    src = dfTmp['{}_src'.format(prop)].unique()

    if val.shape[0] != 1:
        print('No unique {}'.format(prop))

    if src.shape[0] != 1:
        print('No unique {} src'.format(prop))

    if val[0] != '%' and src[0] != '%' and val[0] != '-':
        data = pd.DataFrame([[val[0], src[0]]], columns=['tem', 'larsCode'])

        if prop == 'mlp':
            selectedData.addMeltingPoint(data)
        elif prop == 'blp':
            selectedData.addBoilingPoint(data)
        elif prop == 'tem_cri':
            selectedData.addCriticalTemperature(data)
        elif prop == 'eps':
            selectedData.addPermittivity(data)
        else:
            raise myExceptions.PropertyNotImplemented('selectData::addProperty({})'.format(prop))

    # else:
    #     TODO


def assignCode(dbsConfig):
    """Assigns code to molecules with selected data.
    1. Reads file that maps molecule code to SMILES (if it exists),
    and assigns the codes to molecules.
    2. Create new code for the remaining molecules.

    :param dbsConfig: (dbsConfiguration object) DBS configuration object.
    :return:
    """

    allSelectedData = IO.openSelectedDataFile(dbsConfig)

    fileName = dbsConfig.getOutFileName('codSmilesMap')
    columns = ['code', 'smiles', 'eps']
    codSmilesMap = pd.DataFrame(columns=columns)
    if os.path.exists(fileName):
        codSmilesMap = IO.readCidSmilesMap(fileName, columns)

    # set all codes to empty strings
    for smiles in allSelectedData:
        selectedData = allSelectedData[smiles]
        selectedData.code = ''

    # List of already used molecule codes.
    usedCodes = []

    # assign provided code
    if codSmilesMap.shape[0] != 0:
        # canonicalize SMILES strings.
        for idx, row in codSmilesMap.iterrows():
            smiles = row['smiles']
            smiles = utils.getCanonicalSmiles(smiles)
            codSmilesMap.loc[idx, 'smiles'] = smiles

        # Add code from cod_smiles_map file.
        for smiles in allSelectedData:
            selectedData = allSelectedData[smiles]
            row = codSmilesMap.loc[codSmilesMap['smiles'] == smiles]

            if row.shape[0] != 0:
                code = row['code'].values[0]
                eps = row['eps'].values[0]
                selectedData.code = code

                usedCodes.append(code)
                if np.isnan(eps):
                    print(type(eps), eps)
                if (not np.isnan(eps)) and (not selectedData.eps_src) or (selectedData.eps_src == '-'):
                    selectedData.eps = eps
                    selectedData.eps_src = '-'
                # print('{:6} {} {}'.format(code, smiles, selectedData.eps))

    # add all code in files to usedCodes
    usedCodes.extend(codSmilesMap['code'].values)

    familyCode = dbsConfig.getFamilyCode()
    family = family_utils.getFamily(familyCode)

    # Define code for molecules without one yet.
    for smiles in allSelectedData:
        selectedData = allSelectedData[smiles]
        code = selectedData.code
        # print('{:6} {}'.format(code, smiles, selectedData.eps))

        if not code:
        # if True:
            rdKitMol = Chem.MolFromSmiles(smiles)
            formula = rdMolDescriptors.CalcMolFormula(rdKitMol)
            letter = family.get_letter(smiles)
            # print(formula)

            nC = family.get_num_of_carbons(formula)
            nonC = family.get_num_of_other_atoms_without_H(formula)

            code = utils.createCod(letter, nC, nonC, usedCodes)
            selectedData.code = code

    IO.writeSelectedDataToJson(dbsConfig, allSelectedData)
