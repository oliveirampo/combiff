import matplotlib.pyplot as plt
import numpy as np
import json
import sys
import os


import selectData
import select_points
import myDataStructure


def select_data(fig, ax, x_values, y_valies):
    highlighter = select_points.Highlighter(ax, x_values, y_valies)
    plt.show()
    selected_regions = highlighter.mask

    plt.close(fig)
    return selected_regions


def plotAll(dbsConfig):
    fileName = dbsConfig.getOutFileName('molJsonFile')

    with open(fileName) as jsonFile:
        selectedData = json.load(jsonFile, object_hook=myDataStructure.selectedDataDecoder)

    data = selectedData[0]
    # print(data.mlp, data.blp)
    # print(data)
    print('TODO')


def plotPoint(tem, val, larsCode, marker, color, linewidth):
    plt.plot(tem, val, linewidth=linewidth, marker=marker, c=color, label=larsCode)


def plotStart(mlpVar, mlp, blpVar, blp, tem_cri_var, tem_cri):
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)

    plotRoomTem()
    plotTransitionPoint(mlpVar, mlp)
    plotTransitionPoint(blpVar, blp)
    plotTransitionPoint(tem_cri_var, tem_cri)

    return fig, ax


def plotEnd(fig, prop, cod):
    figName = 'plot/selected_{}_{}.png'.format(prop, cod)
    fig.savefig(figName)
    plt.close(fig)


def plotTransitionPoint(prop, data):
    # plot first value
    values = selectData.getFirstValueOfTransitionPoint(data)

    if len(values) == 2:
        tem = values[0]
        src = values[1]

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

        tem = tem.astype(np.float)
        plt.axvline(x=tem, color=color, linestyle='--', label=label)


# plot hvp at boiling point
def plotHvb(propCod, tables, smiles, dbsEntries, defaultPressure, x_values, y_values, src_values, pre_values,
            fid_values, met_values):

    prop = 'hvb'
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

        marker = dbsRelation.getMarker()
        color = dbsRelation.getColor()

        data = dbsRelation.getData(tab, defaultPressure)

        pre, tem, val, fid, met = dbsRelation.getValues(data)

        plt.plot(tem, val, linewidth=0.0, marker=marker, c=color, label=larsCode)

        for i in range(tem.shape[0]):
            addCode(tem[i], val[i], fid[i])

        saveDataPts(tem, val, larsCode, fid, met, x_values, y_values, src_values, fid_values, met_values)
        for p in pre: pre_values.append(p)


def saveDataPts(X, Y, src, fid, met, x_values, y_values, src_values, fid_values, met_values):
    for i in range(X.shape[0]):
        x_values.append(float(X[i]))
        y_values.append(float(Y[i]))
        fid_values.append(fid[i])
        met_values.append(met[i])

        src_values.append(src)


def markPressure(X, Y, P, defaultPressure):
    for i in range(P.shape[0]):
        x = X[i]
        y = Y[i]
        p = P[i]

        #p = float(p)

        if p > defaultPressure:
            addCode(x, y, '*')


def addCode(x, y, fid):
    plt.annotate(fid, (x, y))


def plotDetails(fig, row, yLabel):
    frm = row['frm']
    nam = row['nam']
    cas = row['cas']
    inchi = row['inchi']

    xLabel = 'T [K]'
    plt.ylabel(yLabel)
    plt.xlabel(xLabel)

    title = '{}\n{} / {}\n{}'.format(frm, nam, cas, inchi)
    plt.title(title)
    plt.xlabel(r'$T\/[K]$')
    plt.legend(numpoints=1, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
    fig.subplots_adjust(right=0.75, top=0.85)


def plotRoomTem():
    plt.axvline(x=298.15, color='black', label='298.15 K')


def insertPic(ax, cod, picDir):
    pngFile = picDir + cod + '.png'
    if os.path.exists(pngFile):
        im = plt.imread(pngFile)
        ax.figure.figimage(im, 460, -60, zorder=-1)


