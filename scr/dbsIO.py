from abc import ABC, abstractmethod
import pandas as pd
import sys
import os


import dbs
import myExceptions


class relationParser(ABC):
    def __init__(self, larsCode, rel, col):
        self.larsCode = larsCode
        self.rel = rel
        self.col = col

    @abstractmethod
    def readRelation(self):
        pass


class parserDefault(relationParser):
    def __init__(self, larsCode, rel, col):
        super(parserDefault, self).__init__(larsCode, rel, col)

    def readRelation(self, path):
        larsCode = self.larsCode
        rel = self.rel
        columns = self.col

        fileName = '{}/{}_{}.rel'.format(path, larsCode, rel)
        if not os.path.exists(fileName):
            raise myExceptions.NoFile(fileName)

        with open(fileName) as dbsFile:
            lines = [row.strip().lower().split() for row in dbsFile.readlines()]

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
                raise myExceptions.WrongNumberofColumns(fileName)

            table.append(row)
            i += 1

    dfTable = pd.DataFrame(table, columns=columns)
    return dfTable


def writeTable(fileName):
    print(fileName)


