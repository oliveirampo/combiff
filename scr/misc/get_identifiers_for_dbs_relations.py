import pubchempy as pcp
import pandas as pd
import time
import sys
import os

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getIdentifiers
import dbsSearch
import utils
import dbs


def search_identifiers(larsCode, dfTable, path):
    if larsCode == 'FR06.6':
        search_identifiers_FR06_6(larsCode, dfTable, path)


def search_identifiers_FR06_6(larsCode, dfTable, path):
    """Iterate though chunks of 1000 rows and search for identifier in PubChem"""

    count = 0

    # add new columns of identifiers
    dfTable['inchi'] = '%'
    dfTable['smiles'] = '%'

    file_name = '{}2/{}_{}.tmp'.format(path, larsCode, 'cpd')
    columns_to_keep = ['cid', 'nam', 'frm', 'cas', 'inchi', 'smiles']
    dfTable.to_csv(file_name, columns=columns_to_keep, sep=' ', index=None)

    for file_count, dfTable in enumerate(pd.read_csv(file_name, sep='\s+', chunksize=1000)):
        file_name = '{}2/{}_{}_{}.rel'.format(path, larsCode, 'cpd', file_count)
        if os.path.exists(file_name):
            continue

        for idx, row in dfTable.iterrows():
            cid = row['cid']
            name = row['nam']
            cas = row['cas']
            formula = row['frm']

            mol_pcp = pcp.get_compounds(name, 'name')

            if len(mol_pcp) != 1:
                print(mol_pcp)
                continue

            mol_pcp = mol_pcp[0]
            mol_cid = mol_pcp.cid
            mol_cas = getIdentifiers.get_cas_pcp(mol_cid)

            mol_name = mol_pcp.iupac_name
            mol_inchi = mol_pcp.inchikey
            mol_formula = mol_pcp.molecular_formula
            mol_smiles = mol_pcp.canonical_smiles

            mol_smiles = utils.getCanonicalSmiles(mol_smiles)

            if cas != '%':
                if cas != mol_cas:
                    sys.exit('{} != {}'.format(cas, mol_cas))

            formula_match = getIdentifiers.match_formula(formula, mol_formula)
            if not formula_match:
                print('{} != {}'.format(formula, mol_formula))
                continue

            if not mol_cas:
                mol_cas = '%'
            if not mol_inchi:
                mol_inchi = '%'
            if not mol_smiles:
                mol_smiles = '%'

            dfTable.loc[idx, 'cas'] = mol_cas
            dfTable.loc[idx, 'inchi'] = mol_inchi
            dfTable.loc[idx, 'smiles'] = mol_smiles

            print(cid, mol_cas, mol_name, mol_smiles)

            count += 1
            if count % 5 == 0:
                time.sleep(1)

        dfTable.to_csv(file_name, sep=' ', index=None)


def main():
    dbsConfigurationFile = '../inp/dbs.conf'
    dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)

    dbsFileName = dbsConfig.getDbsFileName()
    dbsEntries = dbs.getDbsEntries(dbsFileName)

    # path with DBS relations
    path = dbsConfig.getPath()

    for larsCode in dbsEntries:
        factory = dbsEntries[larsCode].getFactory('cpd')
        parser = factory.getParser()
        dfTable = parser.readRelation(path)

        search_identifiers(larsCode, dfTable, path)


if __name__ == '__main__':
    main()
