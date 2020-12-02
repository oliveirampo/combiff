import pandas as pd
import sys


import utils


from molecule import Molecule


def readFieFile(fileName):
    df = pd.read_csv(fileName, sep='\s+', comment='#', names=['nam', 'frm', 'smiles'])
    return df


def readFlsFile(fileName, isomers):
    with open(fileName, 'r') as f:
        lines = f.readlines()

    lines = [row.strip().split() for row in lines]

    molecules = []
    for i in range(len(lines)):
        if not lines[i]:
            continue

        if lines[i][0] == '#':

            # couldn't find a decomposition
            if len(lines[i]) > 2 and lines[i][1] == "couldn\'t":
                smiles = lines[i][6]

                iso = isomers[isomers.smiles == smiles]
                mol = Molecule(smiles, '', '', '', 0, [], 0, '')
                mol.form = iso.frm.iloc[0]
                mol.code = iso.nam.iloc[0]
                mol.enu_code = iso.nam.iloc[0]
                molecules.append(mol)

            else:
                continue

        if lines[i][0][0] == '#':
            continue

        if lines[i][0] == 'NEWMOLECULE':
            pos = i + 1
            mol = readFlsFile_helper(pos, lines)

            smi = mol.smiles
            iso = isomers[isomers.smiles == smi]
            mol.form = iso.frm.iloc[0]
            mol.code = iso.nam.iloc[0]

            molecules.append(mol)

    return molecules


def readFlsFile_helper(pos, data):
    smiles = ''
    name = ''
    code = ''
    cas = ''
    num_frags = 0
    frags = []
    num_links = 0
    links = []

    for i in range(pos, len(data)):
        if data[i][0] == 'ENDMOLECULE':
            mol = Molecule(smiles, name, code, cas, num_frags, frags, num_links, links)
            return mol
        elif data[i][0] == 'NAME':
            smiles = data[i][1]
        elif data[i][0] == 'CODE':
            code = data[i][1]
        elif data[i][0] == 'NUMFRAGS':
            num_frags = data[i][1]
        elif data[i][0] == 'FRAG':
            idx = data[i][1]
            frg = data[i][2]
            frags .append([idx, frg])
        elif data[i][0] == 'NUMLINKS':
            num_links = data[i][1]
        elif data[i][0] == 'LINK':
            frg_1 = data[i][1]
            frg_2 = data[i][2]
            lnk  = data[i][3]
            links.append([frg_1, frg_2, lnk])


def readCiDSmilesFile(fileName):
    df = pd.read_csv(fileName, sep='\s+', header=None, names=['cid', 'smiles'],
        dtype = {'cid': 'Int64', 'smiles': 'str'})

    # canonicalize smiles
    # for idx, row in df.iterrows():
    #     smiles = row['smiles']
    #     smiles = utils.getCanonicalSmiles(smiles)
    #
    #     df.loc[idx, 'smiles'] = str(smiles)

    return df


