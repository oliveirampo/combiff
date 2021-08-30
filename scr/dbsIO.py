"""Module to handle IO of dbs related objets.

Classes:
    relationParser(ABC)
    parserDefault(relationParser)
Methods:
    readTable(i, lines, columns, fileName)
    writeTable(fileName)
"""

from abc import ABC, abstractmethod
import pandas as pd
import os


from scr import myExceptions


class relationParser(ABC):
    """Base class of parser for relation object.
    A relation is a table from a given LARS source and property.

    Attributes:
        larsCode: (str) LARS code.
        rel: (str) Relation code.
        col: (str) String with name/code of every column in relation.
    Functions:
        readRelation()
        removeEmptyValues(var, tab)
    """

    def __init__(self, larsCode, rel, col):
        """Constructs all the necessary attributes for this object.

        :param larsCode: (str) LARS code.
        :param rel: (str) Relation code.
        :param col: (str) String with name/code of every column in relation.
        """

        self.larsCode = larsCode
        self.rel = rel
        self.col = col

    @abstractmethod
    def readRelation(self, path):
        pass

    @staticmethod
    def removeEmptyValues(var, tab):
        """Removes rows of column with name `var` which contain `%` symbol.

        :param var: (str) Column name in table.
        :param tab: (pandas DataFrame) Table with data.
        :return:
            tab: (pandas DataFrame) Table with removed rows.
        """

        if var not in tab.columns:
            return tab

        tab = tab.loc[tab[var] != '%']
        tab = tab.dropna(subset=[var])
        return tab


class parserDefault(relationParser):
    """Default parser for relation object.

    Functions:
        readRelation(path)
        removeEmptyValues(var, tab)
    """

    def readRelation(self, path):
        """Reads relation/table from plain file and returns pandas DataFrame.

        :param path: (str) Path to directory (tmp/) with DBS tables.
        :return:
            dfTable: (pandas DataFrame) Table associated to given relation.
        """

        larsCode = self.larsCode
        rel = self.rel
        columns = self.col

        fileName = '{}/{}_{}.rel'.format(path, larsCode, rel)
        if not os.path.exists(fileName):
            raise myExceptions.NoFile(fileName)

        with open(fileName) as dbsFile:
            lines = [row.strip().split() for row in dbsFile.readlines()]
            # lines = [row.strip().lower().split() for row in dbsFile.readlines()]

        dfTable = pd.DataFrame()

        i = 0
        while i < len(lines):
            row = lines[i]

            if len(row) == 0:
                i += 1

            elif row[0].startswith('#'):
                i += 1

            elif row[0] == '@tab':
                dfTable = readTable(i + 1, lines, columns, fileName)
                i = len(lines)

            else:
                i += 1

        return dfTable


def readTable(i, lines, columns, fileName):
    """Reads file and returns pandas DataFrame

    :param i: (i) Position
    :param lines: (list) Rows from file.
    :param columns: (str) Name of columns of table.
    :param fileName: (str) Name of file.
    :return:
        dfTable: (pandas DataFrame) Table with data.
    """

    columns = columns.strip().split()
    columns = columns[1:]
    table = []

    while i < len(lines):
        row = lines[i]

        if len(row) == 0:
            i += 1

        elif row[0] == '#':
            i = len(lines)

        elif row[0] == '@cap':
            i = len(lines)

        else:
            if len(row) != len(columns):
                print(columns)
                print(row)
                raise myExceptions.WrongNumberofColumns(fileName)

            table.append(row)
            i += 1

    dfTable = pd.DataFrame(table, columns=columns)
    return dfTable
