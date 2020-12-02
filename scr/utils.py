from rdkit import Chem


def getCanonicalSmiles(smiles):
    rdKitMol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(rdKitMol)
    return smiles