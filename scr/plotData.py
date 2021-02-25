"""Module for plotting data.

Methods:
    select_data(fig, ax, x_values, y_valies)
    plotAll(dbsConfig)
    plotPoint(tem, val, larsCode, marker, color, linewidth)
    plotStart(mlpVar, mlp, blpVar, blp, tem_cri_var, tem_cri)
    plotEnd(fig, prop, cod)
    plotTransitionPoint(prop, data)
    plotHvb(propCod, tables, smiles, dbsEntries, defaultPressure, x_values, y_values, src_values, pre_values,
            fid_values, met_values)
    saveDataPts(X, Y, src, fid, met, x_values, y_values, src_values, fid_values, met_values)
    markPressure(X, Y, P, defaultPressure)
    addCode(x, y, fid)
    plotDetails(fig, row, yLabel)
    plotRoomTem()
    insertPic(ax, cod, picDir)
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os

from . import highlight_points
from . import selectData
from . import dbsSearch
from . import dbs
from . import IO


def select_data(fig, ax, x_values, y_values):
    """Selects data points interactively.

    :param fig: (matplotlib.figure.Figure) Figure object.
    :param ax: (matplotlib.axes._subplots.AxesSubplot)
    :param x_values: (numpy.ndarray) X values.
    :param y_values: (numpy.ndarray) Y values.
    :return:
    """

    highlighter = highlight_points.Highlighter(ax, x_values, y_values)
    plt.show()
    selected_regions = highlighter.mask

    plt.close(fig)
    return selected_regions


def plotAll(dbsConfig):
    """Plots all data found + simulated results if provided.

    :param dbsConfig: (dbsConfiguration object) DBS configuration object.
    :return:
    """

    # 00_file.lst
    molListFile = sys.argv[2]
    molList = pd.read_csv(molListFile, sep='\s+', header=None, names=['frm', 'cas', 'nam', 'inchi', 'smiles'])

    allSelectedData = IO.openSelectedDataFile(dbsConfig)

    dbsFileName = dbsConfig.getDbsFileName()
    dbsEntries = dbs.getDbsEntries(dbsFileName)
    defaultPressure = dbsConfig.getDefaultPressure()
    tem_room = dbsConfig.getDefaultTemperature()
    nPoints = dbsConfig.getNumberOfPoints()

    # get tables with data
    properties = dbsConfig.getPropListToBePlotted()
    propCod = dbsConfig.getPropCod()
    tables = dbsSearch.readData(propCod)

    selectData.plotValues(dbsConfig, dbsEntries, defaultPressure, tem_room, nPoints, properties, propCod, tables,
                          molList, allSelectedData, False, False, True, True)


def plotPoint(tem, val, larsCode, marker, color, linewidth):
    """Plots points given x and y coordinates.

    :param tem: (numpy.ndarray) X values.
    :param val: (numpy.ndarray) Y values.
    :param larsCode: (numpy.ndarray) List of LARS codes.
    :param marker: (str) Marker code.
    :param color: (str) Color of points.
    :param linewidth: (float) Line width (should be 0.0).
    :return:
    """

    plt.plot(tem, val, linewidth=linewidth, marker=marker, c=color, label=larsCode)


def plotStart(mlpVar, mlp, blpVar, blp, tem_cri_var, tem_cri):
    """Plots information about molecule:
        - Vertical line for room temperature.
        - Vertical line for melting point.
        - Vertical line for boiling point.
        - Vertical line for critical temperature.

    :param mlpVar: (mlpVar) Melting point code.
    :param mlp: (pandas DataFrame) Table with melting point value.
    :param blpVar: (mlpVar) Boiling point code.
    :param blp: (pandas DataFrame) Table with boiling point value.
    :param tem_cri_var: (mlpVar) Critical temperature code.
    :param tem_cri: (pandas DataFrame) Table with critical temperature.
    :return:
        fig: (matplotlib.figure.Figure) Figure object.
        ax: (matplotlib.axes._subplots.AxesSubplot)
    """

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)

    plotRoomTem()
    plotTransitionPoint(mlpVar, mlp)
    plotTransitionPoint(blpVar, blp)
    plotTransitionPoint(tem_cri_var, tem_cri)

    return fig, ax


def plotEnd(fig, prop, cod):
    """Saves figure into PNG format.

    :param fig: (matplotlib.figure.Figure) Figure object.
    :param prop: (str) Property code.
    :param cod: (cod) Molecule code.
    :return:
    """

    figName = 'plot/selected_{}_{}.png'.format(prop, cod)
    fig.savefig(figName)
    plt.close(fig)


def plotTransitionPoint(prop, data):
    """Plots first value of list for a given transition point.

    :param prop: (str) Property code.
    :param data: (pandas DataFrame) Table with data.
    :return:
    """

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


def plotHvb(propCod, tables, smiles, dbsEntries, defaultPressure, x_values, y_values, src_values, pre_values,
            fid_values, met_values):
    """Plots Vaporization Enthalpy at the boiling point.

    :param propCod: (dict) Dictionary that maps property code to LARDS codes.
    :param tables: (dict) Dictionary of pandas DataFrame that maps property codes to data.
    :param smiles: (str) SMILES string.
    :param dbsEntries: (dict of dbsEntry objects) Dictionary that maps LARS code to dbsEntry objects.
    :param defaultPressure: (float) Default pressure.
    :param x_values: (numpy.ndarray) X values.
    :param y_values: (numpy.ndarray) Y values.
    :param src_values: List of LARS codes.
    :param pre_values: (numpy.ndarray) List of pressure values.
    :param fid_values: (numpy.ndarray) List of annotations.
    :param met_values: (numpy.ndarray) List of methods.
    :return:
    """

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
        dbsRelation = dbsFactory.getRelation()

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
    """Saves values to a unique list.

    :param X: (numpy.ndarray) X values for a given source.
    :param Y: (numpy.ndarray) Y values for a given source.
    :param src: (numpy.ndarray) List of associated LARS codes.
    :param fid: (numpy.ndarray) List of associated annotations.
    :param met: (numpy.ndarray) List of associated methods.
    :param x_values: (numpy.ndarray) List where all X values are saved to.
    :param y_values: (numpy.ndarray) List where all Y values are saved to.
    :param src_values: (numpy.ndarray) List where all LARS codes are saved to.
    :param fid_values: (numpy.ndarray) List where all annotations are saved to.
    :param met_values: (numpy.ndarray) List where all methods are saved to.
    :return:
    """

    for i in range(X.shape[0]):
        x_values.append(float(X[i]))
        y_values.append(float(Y[i]))
        fid_values.append(fid[i])
        met_values.append(met[i])

        src_values.append(src)


def markPressure(X, Y, P, defaultPressure):
    """Add annotation to data point in plot if pressure is higer than default pressure.

    :param X: (numpy.ndarray) X values.
    :param Y: (numpy.ndarray) Y values.
    :param P: (numpy.ndarray) Pressure values.
    :param defaultPressure: (float) Default pressure.
    :return:
    """

    for i in range(P.shape[0]):
        x = X[i]
        y = Y[i]
        p = P[i]

        #p = float(p)

        if p > defaultPressure:
            addCode(x, y, '*')


def addCode(x, y, fid):
    """Adds annotation

    :param x: (float) X coordinate.
    :param y: (float) Y coordinate.
    :param fid: (str) Annotation.
    """

    plt.annotate(fid, (x, y))


def plotDetails(fig, row, yLabel):
    """Plots header of figure.

    :param fig: (matplotlib.figure.Figure) Figure object.
    :param row: (pandas Series) Data to be added to the header of the figure.
    :param yLabel: (str) Y label of plot.
    :return:
    """

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
    """Plots vertical line at room temperature.

    :return:
    """

    plt.axvline(x=298.15, color='black', label='298.15 K')


def insertPic(ax, cod, picDir):
    """Inserts picture to plot.

    :param ax: (matplotlib.axes._subplots.AxesSubplot)
    :param cod: (str) Molecule code.
    :param picDir: (str) Source directory of picture.
    """

    pngFile = picDir + cod + '.png'
    if os.path.exists(pngFile):
        im = plt.imread(pngFile)
        ax.figure.figimage(im, 460, -60, zorder=-1)
