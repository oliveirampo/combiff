from rdkit import Chem


def getCanonicalSmiles(smiles):
    """Returns canonicalized SMILES string by RDKit package.

    :param smiles: (str) SMILES string.
    :return:
        smiles: (str) Canonicalized SMILES string.
    """

    if smiles == '%':
        return smiles

    # print(smiles)
    rdKitMol = Chem.MolFromSmiles(smiles)

    # could not convert smiles to molecule
    if rdKitMol is None:
        print('SMILES not canonicalized: ', smiles)
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


def createCod(letter, nC, nX, usedCodes):
    """Returns molecule code given Family letter code.

    :param letter: (str) One letter family code.
    :param nC: (str) Number of carbon atoms.
    :param nX: (str) Number of other atoms.
    :return:
        code: (str) Molecule code.
    """

    if nC == 10:
        nC = 0

    run = 1
    run_code = '{}{}{}{}'.format(letter, nC, nX, str(run).zfill(2))

    while run_code in usedCodes:
        run += 1
        run_code = '{}{}{}{}'.format(letter, nC, nX, str(run).zfill(2))

    usedCodes.append(run_code)
    return run_code
