from rdkit import Chem


def getCanonicalSmiles(smiles):
    rdKitMol = Chem.MolFromSmiles(smiles)

    # could not convert smiles to molecule
    if rdKitMol == None:
        print(smiles)
        return smiles

    smiles = Chem.MolToSmiles(rdKitMol)
    return smiles


def canonicalizeSmiles(df):
    for idx, row in df.iterrows():
        smiles = row['smiles']
        smiles = getCanonicalSmiles(smiles)
        df.loc[idx, 'smiles'] = smiles

    return df


def createCod(cur_family, frm, letter, run, used):
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