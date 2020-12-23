import matplotlib.pyplot as plt
import json
import sys
import os


import selectData
import select_points
import myDataStructure


def select_data(fig, ax, x_values, y_valies):
    highlighter = select_points.Highlighter(ax, x_values, y_valies)
    # plt.show()
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
    fig = plt.figure()
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

        plt.axvline(x=tem, color=color, linestyle='--', label=label)


# plot hvp at boiling point
def plotHvb(propCod, tables, smiles, dbsEntries, defaultPressure, x_values, y_values, src_values, pre_values):
    prop = 'hvb'

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


def saveDataPts(X, Y, src, x_values, y_values, src_values):
    for x in X:
        x_values.append(float(x))
        src_values.append(src)
    for y in Y:
        y_values.append(float(y))


def markPressure(X, Y, P, defaultPressure):
    for i in range(P.shape[0]):
        x = X[i]
        y = Y[i]
        p = P[i]

        if p > defaultPressure:
            addCode(x, y, '*')


def addCode(x, y, fid):
    plt.annotate(fid, (x, y))


def getYLabel(prop):
    labels = {
        'dns': r'$\rho_{liq} \/ [kg \cdot m^{-1}]$',
        'hvp': r'$\Delta H_{vap} \/ [kJ \cdot mol^{-1}]$'
    }
    yLab = labels[prop]
    return yLab


def plotDetails(fig, row, prop):
    cod = row['cod']
    frm = row['frm']
    nam = row['nam']
    cas = row['cas']
    inchi = row['inchi']

    yLabel = getYLabel(prop)
    xLabel = 'T [K]'
    plt.ylabel(yLabel)
    plt.xlabel(xLabel)

    title = '{} / {}\n{} / {}\n{}'.format(cod, frm, nam, cas, inchi)
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


