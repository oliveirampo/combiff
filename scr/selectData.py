import pandas as pd
import numpy as np
import json
import sys


import myDataStructure
import myExceptions
import dbsSearch
import plotData
import dbs
import IO


def manualSelection(dbsConfig):
    nArgs = len(sys.argv)
    if nArgs != 3:
        raise myExceptions.ArgError(nArgs, 3)

    # 00_file.lst
    molListFile = sys.argv[2]
    molList = pd.read_csv(molListFile, sep='\s+', header=None, names=['cod', 'frm', 'cas', 'nam', 'inchi', 'smiles'])

    dbsFileName = dbsConfig.getDbsFileName()
    dbsEntries = dbs.getDbsEntries(dbsFileName)

    # get tables with data
    properties = dbsConfig.getProps()
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
            IO.writeToJson(dbsConfig, allSelectedData)


def getTransitionPoint(propName, defaultPressure, dbsConfig, dbsEntries, propCod, tables, smiles):
    prop = dbsConfig.getVariable(propName)

    data = []
    for larsCode in propCod[prop]:
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




