from rdkit import Chem


def getCanonicalSmiles(smiles):
    """Returns canonicalized SMILES string by RDKit package.

    :param smiles: (str) SMILES string.
    :return:
        smiles: (str) Canonicalized SMILES string.
    """

    rdKitMol = Chem.MolFromSmiles(smiles)

    # could not convert smiles to molecule
    if rdKitMol is None:
        print(smiles)
        return smiles

    smiles = Chem.MolToSmiles(rdKitMol)
    return smiles


def canonicalizeSmiles(df):
    """Returns table with canonizalized SMILES strings.

    :param df: (pandas DataFrame) Table.
    :return:
        df: (pandas DataFrame) Table with canonicalized SMILES strings..
    """

    for idx, row in df.iterrows():
        smiles = row['smiles']
        smiles = getCanonicalSmiles(smiles)
        df.loc[idx, 'smiles'] = smiles

    return df


def createCod(cur_family, frm, letter, run, used):
    """Returns molecule code given Family letter code.

    :param cur_family:
    :param frm:
    :param letter:
    :param run:
    :param used:
    :return:
    """
    # TODO

    nX = cur_family.get_num_of_other_atoms(frm)
    nC = cur_family.get_num_of_carbons(frm)

    if nC == 10:
        nC = 0

    run_code = '{}{}{}{}'.format(letter, nC, nX, str(run).zfill(2))

    while run_code in used:
        run += 1
        run_code = '{}{}{}{}'.format(letter, nC, nX, str(run).zfill(2))

    used.append(run_code)
    return run_code
