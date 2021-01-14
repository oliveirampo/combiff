import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os


import myDataStructure
import myExceptions
import dbsSearch
import plotData
import utils
import dbs
import IO


def manualSelection(dbsConfig):
    nArgs = len(sys.argv)
    if nArgs != 3:
        raise myExceptions.ArgError(nArgs, 3)

    # 00_file.lst
    molListFile = sys.argv[2]
    molList = pd.read_csv(molListFile, sep='\s+', header=None, names=['frm', 'cas', 'nam', 'inchi', 'smiles'])

    dbsFileName = dbsConfig.getDbsFileName()
    dbsEntries = dbs.getDbsEntries(dbsFileName)

    # get tables with data
    properties = dbsConfig.getPropList()
    propCod = dbsConfig.getPropCod()
    tables = dbsSearch.readData(propCod)

    defaultPressure = dbsConfig.getDefaultPressure()
    tem_room = dbsConfig.getDefaultTemperature()
    nPoints = dbsConfig.getNumberOfPoints()

    allSelectedData = []

    for prop in properties:
        for idx, row in molList.iterrows():
            smiles = row['smiles']
            selectedData = myDataStructure.SelectedData(smiles)

            mlpVar, mlp = getTransitionPoint('meltingPoint', defaultPressure, dbsConfig, dbsEntries, propCod, tables, smiles)
            blpVar, blp = getTransitionPoint('boilingPoint', defaultPressure, dbsConfig, dbsEntries, propCod, tables, smiles)
            tem_cri_var, tem_cri = getTransitionPoint('criticalTemperature', defaultPressure, dbsConfig, dbsEntries, propCod, tables, smiles)

            selectedData.addMeltingPoint(mlp)
            selectedData.addBoilingPoint(blp)
            selectedData.addCriticalTemperature(tem_cri)

            fig, ax = plotData.plotStart(mlpVar, mlp, blpVar, blp, tem_cri_var, tem_cri)

            x_values = []
            y_values = []
            pre_values = []
            src_values = []

            if dbsConfig.isHvp(prop):
                plotData.plotHvb(propCod, tables, smiles, dbsEntries, defaultPressure, x_values, y_values, src_values, pre_values)

            for larsCode in propCod[prop]:
                if larsCode not in tables[prop]:
                    continue

                tab = tables[prop][larsCode]
                tab = tab.loc[tab['smiles'] == smiles]
                if tab.shape[0] == 0:
                    continue

                dbsEntry = dbsEntries[larsCode]
                dbsFactory = dbsEntry.getFactory(prop)
                dbsRelation = dbsFactory.createRelation()
                data = dbsRelation.getData(tab, defaultPressure)

                pre = data[:, 0]
                tem = data[:, 1]
                val = data[:, 2]

                marker = dbsRelation.marker
                color = dbsRelation.color
                plotData.plotPoint(tem, val, larsCode, marker, color, linewidth=0.0)

                plotData.saveDataPts(tem, val, larsCode, x_values, y_values, src_values)
                for p in pre: pre_values.append(p)

                equation = dbsFactory.createEquation()
                X, Y, note = equation.getData(tem_room, nPoints, tab, dbsRelation.tem_convert)
                if len(X) > 0:
                    plotData.plotPoint(X, Y, larsCode, marker, color, linewidth=1.0)
                    plotData.saveDataPts(X, Y, larsCode, x_values, y_values, src_values)
                    plotData.addCode(X[0], Y[0], note[0])

                    for x in X: pre_values.append(defaultPressure)

            x_values, y_values, pre_values, src_values = checkArrays(prop, x_values, y_values, pre_values, src_values)

            plotData.plotDetails(fig, row, prop)

            addVaporPressure(propCod, blp, tables, smiles, dbsEntries, x_values, pre_values, src_values)

            plotData.markPressure(x_values, y_values, pre_values, defaultPressure)

            selectedPointsMask = []
            if len(x_values) != 0: selectedPointsMask = plotData.select_data(fig, ax, x_values, y_values)

            selectedData.addProperty(prop, selectedPointsMask, x_values, y_values, pre_values, src_values)
            allSelectedData.append(selectedData)
            # IO.writeSelectedDataToJson(dbsConfig, allSelectedData)

            plt.close(fig)



def getTransitionPoint(propName, defaultPressure, dbsConfig, dbsEntries, propCod, tables, smiles):
    prop = dbsConfig.getVariable(propName)

    data = []
    for larsCode in propCod[prop]:
        if larsCode not in tables[prop]:
            continue

        tab = tables[prop][larsCode]
        tab = tab.loc[tab['smiles'] == smiles]

        if tab.shape[0] == 0:
            continue

        dbsEntry = dbsEntries[larsCode]
        dbsFactory = dbsEntry.getFactory(prop)
        dbsRelation = dbsFactory.createRelation()
        dat = dbsRelation.getData(tab, defaultPressure)

        values = dat.shape[0] * [[0.0, larsCode]]

        for i in range(dat.shape[0]):
            values[i][0] = dat[i, 0]
        data.extend(values)

    data = pd.DataFrame(data, columns=['tem', 'larsCode'])
    return prop, data


def addVaporPressure(propCod, blp, tables, smiles, dbsEntries, x_values, pre_values, src_values):
    prop = 'pvp'
    if not prop in propCod:
        return

    maxBlp = getHighestTransitionPoint('tem', blp)

    for larsCode in propCod[prop]:
        tab = tables[prop][larsCode]
        tab = tab.loc[tab['smiles'] == smiles]
        if tab.shape[0] == 0:
            continue

        dbsEntry = dbsEntries[larsCode]
        dbsFactory = dbsEntry.getFactory(prop)

        equation = dbsFactory.createEquation()

        for i in range(len(src_values)):
            src = src_values[i]

            if src == 'YA14.6':
                tem = x_values[i]

                if tem > maxBlp:
                    val = equation.getDataAt(tem, tab)
                    if len(val) == 1:
                        val = val[0]
                        pre_values[i] = val


def getFirstValueOfTransitionPoint(data):
    # plot first value
    if data.shape[0] != 0:
        row = data.iloc[0]
        tem = row['tem']
        src = row['larsCode']
        return [tem, src]
    return []


def getHighestTransitionPoint(var, data):
    maxVal = data[var].max()
    maxVal = float(maxVal)
    return maxVal


def checkArrays(prop, x_values, y_values, pre_values, src_values):
    x_values = np.asarray(x_values)
    y_values = np.asarray(y_values)
    pre_values = np.asarray(pre_values)
    src_values = np.asarray(src_values)

    if len(x_values) != len(y_values):
        print('len(tem) != len({})\n'.format(prop))
        sys.exit(1)
    if len(x_values) != len(pre_values):
        print('len(tem) != len(pre) for prop = {}\n'.format(prop))
        sys.exit(1)

    return x_values, y_values, pre_values, src_values


# get selected data from old format file
# and add to already available output file
def getSelectedDataForProperty(dbsConfig):
    fieFile = sys.argv[3]
    isomers = IO.readFieFile(fieFile)

    oldMolListFileName = sys.argv[4]    # 01_mol.lst

    molList = pd.read_csv(oldMolListFileName, sep='\s+',
              names=['cod', 'frm', 'cas', 'nam', 'inchi', 'smiles', 'old_cod', 'x', 'y', 'z'],
              usecols=['cod', 'inchi', 'smiles'])

    propCod = dbsConfig.getPropCod()
    allSelectedData = IO.openSelectedDataFile(dbsConfig)

    for prop in propCod:
        print(prop)
        if prop == 'diffus':
            oldSelectedDataFileName = 'old_files/{}.src'.format('self_dif')
        elif prop == 'etd':
            oldSelectedDataFileName = 'old_files/{}.src'.format('vis_d')
        else:
            oldSelectedDataFileName = 'old_files/{}.src'.format(prop)

        if not os.path.exists(oldSelectedDataFileName):
            continue

        srcData = pd.read_csv(oldSelectedDataFileName, sep='\s+', names=['cod', 'pre', 'tem', 'val', 'src'])
        srcData = srcData.loc[srcData['val'].astype(str) != '%']

        mergedData = pd.merge(srcData, molList, how='inner', on='cod')

        # get smiles in common
        isomers = utils.canonicalizeSmiles(isomers)
        mergedData = utils.canonicalizeSmiles(mergedData)
        mergedData = pd.merge(isomers, mergedData, how='inner', on='smiles')

        for selectedData in allSelectedData:
            smiles = selectedData.smiles

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

            if not selectedData.hasProperty(prop):
                property = myDataStructure.Property(prop, pre, tem, val, src)
                selectedData.addPropertyHelper(prop, property)

    IO.writeSelectedDataToJson(dbsConfig, allSelectedData)


# get selected data from old format file
# for dns and hvp
def getSelectedData(dbsConfig):
    fieFile = sys.argv[3]
    isomers = IO.readFieFile(fieFile)

    oldMolListFileName = sys.argv[4]    # 01_mol.lst
    oldSelectedDataFileName = 'old_files/mol.src'

    molList = pd.read_csv(oldMolListFileName, sep='\s+',
        names=['cod', 'frm', 'cas', 'nam', 'inchi', 'smiles', 'old_cod', 'x', 'y', 'z'],
        usecols=['cod', 'inchi', 'smiles'])

    srcData = pd.read_csv(oldSelectedDataFileName, sep='\s+',
        names = ['cod', 'nam', 'cas', 'pre', 'tem', 'dns', 'dns_src', 'hvp', 'hvp_src',
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
    allSelectedData = []
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
            pre = data['pre'].values
            tem = data['tem'].values
            dns = data['dns'].values
            src = data['dns_src'].values
            appendData = True

            property = myDataStructure.Property('dns', pre, tem, dns, src)
            selectedData.addPropertyHelper('dns', property)

        data = dfTmp[dfTmp['hvp'].astype(str) != '%']
        data = data[data['hvp'].astype(str) != '-']
        if data.shape[0] != 0:
            pre = data['pre'].values
            tem = data['tem'].values
            hvp = data['hvp'].values
            src = data['hvp_src'].values
            appendData = True

            property = myDataStructure.Property('hvp', pre, tem, hvp, src)
            selectedData.addPropertyHelper('hvp', property)

        if appendData:
            allSelectedData.append(selectedData)
            IO.writeSelectedDataToJson(dbsConfig, allSelectedData)


def addProperty(prop, dfTmp, selectedData):
    val = dfTmp[prop].unique()
    src = dfTmp['{}_src'.format(prop)].unique()

    if val.shape[0] != 1:
        print('No unique {}'.format(prop))

    if src.shape[0] != 1:
        print('No unique {} src'.format(prop))

    if val[0] != '%' and src[0] != '%' and val[0] != '-' and src[0] != '-':
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
