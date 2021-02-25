"""Writes identifiers to txt file.

Methods:
    run()
    selectMolecules(isomer, enuData)
"""

import json
import sys
import os


from .molecule import moleculeDecoder
from . import myExceptions
from . import utils
from . import IO


def run():
    """Writes identifiers to txt file."""

    nArgs = len(sys.argv)
    if nArgs != 3 and nArgs != 4:
        raise myExceptions.ArgError('3 or 4', nArgs)

    name = sys.argv[2]
    name = name.split('/')[-1]

    enuMolFile = 'out/{}.json'.format(name)
    if not os.path.exists(enuMolFile):
        raise myExceptions.NoFile(enuMolFile)

    with open(enuMolFile) as jsonFile:
        enuData = json.load(jsonFile, object_hook=moleculeDecoder)

    if len(sys.argv) == 4:
        fieFile = sys.argv[3]
        isomers = IO.readFieFile(fieFile)
        isomers = utils.canonicalizeSmiles(isomers)

        enuData = selectMolecules(isomers, enuData)

    outFileName = 'out/00_{}.lst'.format(name)
    out = open(outFileName, 'w')

    for smiles in enuData:
        mol = enuData[smiles]
        # frm = mol.form
        frm = mol.form_pcp

        cas = mol.cas
        smiles = mol.smiles

        if not cas:
            cas = '%'

        name = mol.name_pcp
        if not name:
            name = '%'

        inchi = mol.inchi_pcp

        name = name.replace(' ', '_')
        out.write('{:10} {:12} {:40} {:30} {}\n'
                  .format(frm, cas, name, inchi, smiles))

    out.close()


def selectMolecules(isomer, enuData):
    """Returns data of only molecules specified in isomer DataFrame.

    :param isomer: (pandas DataFrame) List of constitutional isomers together with molecular formula and SMILES.
    :param enuData: (dict) Molecule data (identifiers).
    :return:
        data: (dict) Selected molecule data (identifiers).
    """

    data = {}
    for idx, row in isomer.iterrows():
        smiles = row['smiles']

        if smiles in enuData:
            mol = enuData[smiles]
            data[smiles] = mol

    return data
