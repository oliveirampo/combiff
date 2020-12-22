import matplotlib.pyplot as plt
import pandas as pd
import configparser
import numpy as np
import sys
import os


import myExceptions
import dbsSearch
import dbs


def plotAll(dbsConfig):
    nArgs = len(sys.argv)
    if nArgs != 3:
        raise myExceptions.ArgError(nArgs, 3)

    # 00_file.lst
    molListFile = sys.argv[2]
    molList = pd.read_csv(molListFile, sep='\s+', header=None, names=['cod', 'frm', 'cas', 'nam', 'inchi', 'smiles'])
    # print(molList)

    dbsFileName = dbsConfig.getDbsFileName()
    dbsEntries = dbs.getDbsEntries(dbsFileName)

    # get tables with data
    properties = dbsConfig.getProps()
    propCod = dbsConfig.getPropCod()
    tables = dbsSearch.readData(propCod)

    defaultPressure = dbsConfig.getDefaultPressure()
    tem_room = dbsConfig.getDefaultTemperature()
    nPoints = dbsConfig.getNumberOfPoints()

    for prop in properties:
        for idx, row in molList.iterrows():
            cod = row['cod']
            frm = row['frm']
            nam = row['nam']
            cas = row['cas']
            inchi = row['inchi']
            smiles = row['smiles']

            fig, ax = plotStart()
            plotRoomTem()
            plotTransitionPoint('meltingPoint', defaultPressure, dbsConfig, dbsEntries, propCod, tables, smiles)
            plotTransitionPoint('boilingPoint', defaultPressure, dbsConfig, dbsEntries, propCod, tables, smiles)
            plotTransitionPoint('criticalTemperature', defaultPressure, dbsConfig, dbsEntries, propCod, tables, smiles)

            x_values = []
            y_values = []
            src_values = []
            pre_values = []

            if dbsConfig.isHvp(prop):
                plotHvb('hvb', propCod, tables, smiles, dbsEntries, defaultPressure, x_values, y_values, src_values, pre_values)

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
                plt.plot(tem, val, linewidth=0.0, marker=marker, c=color, label=larsCode)
                saveDataPts(tem, val, larsCode, x_values, y_values, src_values)
                for p in pre: pre_values.append(p)

                equation = dbsFactory.createEquation()
                X, Y, note = equation.getData(tem_room, nPoints, tab, dbsRelation.tem_convert)
                if len(X) > 0:
                    plt.plot(X, Y, linewidth=1.0, marker=marker, color=color, label=larsCode)
                    saveDataPts(X, Y, larsCode, x_values, y_values, src_values)
                    addCode(X[0], Y[0], note[0])

                    # TODO - add pvp
                    for x in X: pre_values.append(defaultPressure)

            # TODO - add pvp

            yLabel = getYLabel(prop)
            xLabel = 'T [K]'
            plt.ylabel(yLabel)
            plt.xlabel(xLabel)

            plotDetails(fig, cod, frm, nam, cas, inchi)

            x_values = np.asarray(x_values)
            y_values = np.asarray(y_values)
            src_values = np.asarray(src_values)
            pre_values = np.asarray(pre_values)

            if len(x_values) != len(y_values):
                print('len(tem) != len({})\n'.format(prop))
                sys.exit(1)
            if len(x_values) != len(pre_values):
                print('len(tem) != len(pre) for prop = {}\n'.format(prop))
                sys.exit(1)

            # plt.show()
            figName = 'plot/all_{}_{}.png'.format(prop, cod)
            fig.savefig(figName)
            plt.close(fig)


# plot hvp at boiling point
def plotHvb(prop, propCod, tables, smiles, dbsEntries, defaultPressure, x_values, y_values, src_values, pre_values):
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
        plt.plot(tem, val, linewidth=0.0, marker=marker, c=color, label=larsCode)

        saveDataPts(tem, val, larsCode, x_values, y_values, src_values)
        for p in pre: pre_values.append(p)


def addCode(x, y, fid):
    plt.annotate(fid, (x, y))


def saveDataPts(X, Y, src, x_values, y_values, src_values):
    for x in X:
        x_values.append(float(x))
        src_values.append(src)
    for y in Y:
        y_values.append(float(y))


def getYLabel(prop):
    labels = {
        'dns': r'$\rho_{liq} \/ [kg \cdot m^{-1}]$',
        'hvp': r'$\Delta H_{vap} \/ [kJ \cdot mol^{-1}]$'
    }
    yLab = labels[prop]
    return yLab


def plotStart():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    return fig, ax


def plotRoomTem():
    plt.axvline(x=298.15, color='black', label='298.15 K')


def plotTransitionPoint(propName, defaultPressure, dbsConfig, dbsEntries, propCod, tables, smiles):
    try:
        prop = dbsConfig.getVariable(propName)
    except configparser.NoOptionError:
        raise myExceptions.VariableNotDefined(propName)

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

    # plot first value
    if data.shape[0] != 0:
        row = data.iloc[0]
        tem = row['tem']
        src = row['larsCode']

        if prop == 'mlp':
            color = 'blue'
            label = 'Tm (' + src + ')'
        elif prop == 'blp':
            color = 'lime'
            label = 'Tb (' + src + ')'
        elif prop == 'tem_cri':
            color = 'red'
            label = 'Tc (' + src + ')'
        else:
            print('Variable not defined: {}'.format(prop))
            sys.exit(123)

        plt.axvline(x=tem, color=color, linestyle='--', label=label)


def plotDetails(fig, cod, frm, nam, cas, stdInchiKey):
    title = '{} / {}\n{} / {}\n{}'.format(cod, frm, nam, cas, stdInchiKey)
    plt.title(title)
    plt.xlabel(r'$T\/[K]$')
    plt.legend(numpoints=1, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
    fig.subplots_adjust(right=0.75, top=0.85)


def insertPic(ax, cod, picDir):
    pngFile = picDir + cod + '.png'
    if os.path.exists(pngFile):
        im = plt.imread(pngFile)
        ax.figure.figimage(im, 460, -60, zorder=-1)


