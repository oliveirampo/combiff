import json
import sys
import os


import myExceptions
import family_utils
from molecule import moleculeDecoder


def run():
    nArgs = len(sys.argv)
    if nArgs != 4:
        raise myExceptions.ArgError(nArgs, 4)

    name = sys.argv[2]
    name = name.split('/')[-1]
    familyCod = sys.argv[3]

    enuMolFile = 'inp/{}.json'.format(name)
    if not os.path.exists(enuMolFile):
        raise myExceptions.NoFile(enuMolFile)

    with open(enuMolFile) as jsonFile:
        enuData = json.load(jsonFile, object_hook=moleculeDecoder)

    outFileName = 'inp/00_{}.lst'.format(name)
    out = open(outFileName, 'w')
    used = []

    cur_family = family_utils.getFamily(familyCod)
    letter = cur_family.letter

    for smiles in enuData:
        run = 0
        mol = enuData[smiles]
        frm = mol.form
        cas = mol.cas
        smiles = mol.smiles

        nX = cur_family.get_num_of_other_atoms(frm)
        nC = cur_family.get_num_of_carbons(frm)

        if nC == 10:
            nC = 0

        run_code = '{}{}{}{}'.format(letter, nC, nX, str(run).zfill(2))

        while run_code in used:
            run += 1
            run_code = '{}{}{}{}'.format(letter, nC, nX, str(run).zfill(2))

        used.append(run_code)

        if not cas:
            cas = '%'

        name = mol.name_pcp
        if not name:
            name = '%'

        inchi = mol.inchi_pcp

        name = name.replace(' ', '_')
        out.write('{:5} {:10} {:10} {:40} {:30} {}\n'
            .format(run_code, frm, cas, name, inchi, smiles))

    out.close()