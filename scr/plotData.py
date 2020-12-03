import matplotlib.pyplot as plt
import pandas as pd
import sys


import myExceptions
import dbsSearch
import dbs


def plotAll():
    nArgs = len(sys.argv)
    if nArgs != 3:
        raise myExceptions.ArgError(nArgs, 3)

    # 00_file.lst
    molListFile = sys.argv[2]
    molList = pd.read_csv(molListFile, sep='\s+', header=None, names=['cod', 'frm', 'cas', 'nam', 'inchi', 'smiles'])
    # print(molList)

    dbsFileName = 'inp/dbs.json'
    dbsEntries = dbs.getDbsEntries(dbsFileName)

    # get tables with data
    tables = dbsSearch.readData()

    propCod = dbs.dbsEntry.propCod
    for prop in propCod:

        for idx, row in molList.iterrows():
            # cas = row['cas']
            # nam = row['nam']
            cod = row['cod']
            smiles = row['smiles']

            fig, ax = plotStart()
            plotRoomTem()
            # TODO - plotBlp, plotMlp, plotTemCri

            for larsCode in propCod[prop]:
                tab = tables[prop][larsCode]
                tab = tab.loc[tab['smiles'] == smiles]
                if tab.shape[0] == 0:
                    continue

                dbsEntry = dbsEntries[larsCode]
                dbsFactory = dbsEntry.getFactory(prop)
                dbsRelation = dbsFactory.createRelation()
                data = dbsRelation.getData(tab)

                pre = data[:, 0]
                tem = data[:, 1]
                val = data[:, 2]
                marker = dbsRelation.marker
                color = dbsRelation.color
                pts = plt.plot(tem, val, linewidth=0.0, marker=marker, c=color, label=larsCode)


            yLabel = getYLabel(prop)
            xLabel = 'T [K]'
            plt.ylabel(yLabel)
            plt.xlabel(xLabel)

            # plt.show()
            figName = 'plot/all_{}_{}.png'.format(prop, cod)
            fig.savefig(figName)
            plt.close(fig)


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


