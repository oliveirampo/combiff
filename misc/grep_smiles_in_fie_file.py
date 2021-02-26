import pandas as pd
import numpy as np
import sys
import os

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from scr import IO
from scr import utils




def main():
    """Greps the SMILES strings listed in molList file from the fie file.

    :return:
    """

    molListFile = sys.argv[1]
    fieFile = sys.argv[2]
    outFile = sys.argv[3]

    isomers = IO.readFieFile(fieFile)
    isomers = utils.canonicalizeSmiles(isomers)

    molList = pd.read_csv(molListFile, sep='\s+',
                         names=['code', 'formula', 'cas', 'name', 'inchi', 'smiles', 'K', 'X', 'Y', 'Z'])
    molList = utils.canonicalizeSmiles(molList)

    df = pd.merge(molList, isomers, on='smiles')

    columns = ['nam', 'frm', 'smiles']
    if df.shape[0] != molList.shape[0]:
        print('ERROR: Some SMILES were not matched.')
        # sys.exit('ERROR: Some SMILES were not matched.')

        df = df[columns]
        # show duplicates
        duplicated_smiles = np.unique(df[df.duplicated(['smiles'], keep=False)]['smiles'])
        for smiles in duplicated_smiles:
            print(smiles)

    df.to_csv(outFile, columns=columns, index=False, header=False, sep='\t')


if __name__ == '__main__':
    main()