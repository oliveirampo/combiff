"""Gets identifiers (name, cas, standard inchikey) from pubchem given SMILES strings.
    - If a file with CID and smiles from pubchem is given,
    then the CID is used directly to retrieve the molecule from pubchem.
    - If an empty file is given,
    then the smiles is used to search for the molecule in the database.
In both cases RDKit is used to canonicalize the smiles strings.

Methods:
"""

from pubchempy import BadRequestError
from collections import OrderedDict
import pubchempy as pcp
import time
import json
import sys
import os
import re

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from molecule import moleculeEncoder
from molecule import moleculeDecoder
from molecule import Molecule
import myExceptions
import utils
import IO


def run(dbsConfig):
    """Gets identifiers (name, cas, standard inchikey) from pubchem given SMILES strings.
    1. Read input files.
        - fieFile (list of constitutional isomers along with molecular formula ans SMILES strings.
        - empty file or file that maps CID from PubChem to SMILES strings.
    2. Canonicalize SMILES using RDKit.
    3. Match SMILES from fieFile to SMILES from PubChem and get CID from PubChem.
        - If a file with CID and smiles from pubchem is given,
        then the CID is used directly to retrieve the molecule from pubchem.
        - If an empty file is given,
        then the smiles is used to search for the molecule in the database.
    4. Use CID from PubChem to download molecule and extract identifiers.
        - Check if formula match.

    :param dbsConfig: (dbsConfiguration object) DBS configuration object.
    :return:
    """

    nArgs = len(sys.argv)
    if nArgs < 4:
        raise myExceptions.ArgError(4, nArgs)

    # flsFile = sys.argv[2]
    fieFile = sys.argv[2]
    cidSmilesFile = sys.argv[3]

    isomers = IO.readFieFile(fieFile)
    # molecules = IO.readFlsFile(flsFile, isomers)
    molecules = getMolecules(isomers)

    print('Canonicalizing SMILES')
    canonicalizeSmiles(molecules)
    # this will take a while
    cidSmiles = IO.readCiDSmilesFile(cidSmilesFile)

    data = OrderedDict()
    # check if json/pickle file exists
    # molDataFile = 'out/enuMolData.json'
    molDataFile = dbsConfig.getOutFileName('enuDataFile')
    if os.path.exists(molDataFile):
        with open(molDataFile) as jsonFile:
            data = json.load(jsonFile, object_hook=moleculeDecoder)

    print('Downloading Molecules')
    if cidSmiles.shape[0] != 0:
        matchMol(cidSmiles, molecules, data, molDataFile)
    else:
        downloadMol(molecules, data, molDataFile)


def getMolecules(isomers):
    """Returns list of molecules created from table with SMILES.

    :param isomers: (pandas DataFrame) List of constitutional isomers together with molecular formula and SMILES.
    :return:
        molecules: (list) List of Molecules.
    """

    molecules = []
    for idx, row in isomers.iterrows():
        smiles = row['smiles']
        form = row['frm']
        nam = row['nam']

        mol = Molecule(smiles, '', nam, '', 0, [], 0, '')
        mol.form = form
        molecules.append(mol)

    return molecules


def canonicalizeSmiles(molecules):
    """ Canonicalize SMILES strings from list of molecules.

    :param molecules: (list) List of molecules.
    """

    for mol in molecules:
        smiles = mol.smiles
        smiles = utils.getCanonicalSmiles(smiles)
        mol.smiles = smiles


def matchMol(cidSmiles, molecules, data, molDataFile):
    """Extracts identifier of molecules.
    1. Match smiles from fieFile to smiles from file that maps CID from PubChem to SMILES strings.
    2. Use the CID to download molecule from PubChem using the pubchempy library.
    3. Extract identifiers.
    4. Check consistency between formulas.

    :param cidSmiles: (pandas DataFrame) Table that maps CID from PubChem to SMILES strings.
    :param molecules: (list) Molecules.
    :param data: (OrderedDict) Dictionary of downloaded data. Additional data is appended to it.
    :param molDataFile: Output file to it data will be written.
    :return:

    There is higher chance of matching identifiers then using the downloadMol() method
    because the smiles to be matched are canonicalized with RDKit first.

    Wait 1 sec after searching for 5 molecules
    """

    count = 0
    for mol in molecules:
        code = mol.enu_code
        smiles = mol.smiles

        if smiles in data:
            print('\t{} ({}) already included'.format(code, smiles))
            continue

        print(mol)

        cid_pcp = cidSmiles[cidSmiles['smiles'] == smiles]['cid']

        if cid_pcp.shape[0] == 1:
            cid_pcp = cid_pcp.values[0]
            cid_pcp = str(cid_pcp)
            mol_pcp = pcp.get_compounds(cid_pcp, 'cid')

            if len(mol_pcp) == 0:
                continue
            if len(mol_pcp) != 1:
                print('More than 1 match for: {} {}'.format(mol_pcp, smiles))
                sys.exit(1)

            mol_pcp = mol_pcp[0]
            if mol_pcp.cid is None:
                continue

            getPcpData(mol_pcp, smiles, mol, data, molDataFile)

            count += 1
            if count % 5 == 0:
                time.sleep(1)


def downloadMol(molecules, data, molDataFile):
    """Extracts identifier of molecules.
    1. Use smiles from fieFile to search for molecules in PubChem using pubchempy library.
    2. Extract identifiers.
    3. Check consistency between formulas.

    :param molecules: (list) Molecules.
    :param data: (OrderedDict) Dictionary of downloaded data. Additional data is appended to it.
    :param molDataFile: Output file to it data will be written.

    Wait 1 sec after searching for 5 molecules
    """

    count = 0

    for mol in molecules:
        code = mol.enu_code
        smiles = mol.smiles

        # canonicalize smiles with RDKit
        smiles = utils.getCanonicalSmiles(smiles)
        mol.smiles = smiles

        if smiles in data:
            print('\t{} ({}) already included'.format(code, smiles))
            continue

        print(mol)

        try:
            mol_pcp = pcp.get_compounds(smiles, 'smiles')
        except BadRequestError:
            print('\n\tBadRequesteError for SMILES = {}\n'.format(smiles))
            continue

        if len(mol_pcp) == 0:
            continue
        if len(mol_pcp) != 1:
            print('More than 1 match for: {} {}'.format(mol_pcp, smiles))
            sys.exit(1)

        mol_pcp = mol_pcp[0]
        if mol_pcp.cid is None:
            continue

        getPcpData(mol_pcp, smiles, mol, data, molDataFile)

        count += 1
        if count % 5 == 0:
            time.sleep(1)


def getPcpData(mol_pcp, smiles, mol, data, molDataFile):
    """Extracts information from pubchempy.Compound and save it to json file.

    :param mol_pcp: (pubchempy.Compound) Molecule from PubChem.
    :param smiles: (str) SMILES string.
    :param mol: (molecule.Molecule) Molecule from fieFile.
    :param data: (OrderedDict) Dictionary of downloaded data. Additional data is appended to it.
    :param molDataFile: Output file to it data will be written.
    :return:
    """

    name = mol_pcp.iupac_name

    if not name:
        name = getFirstSynonym(mol_pcp)
    add_other_names(mol, name)

    match = compare(mol, mol_pcp)
    if not match:
        return
    add_other_names(mol, mol_pcp)

    data[smiles] = mol
    # save data to json file
    with open(molDataFile, 'w') as out:
        json.dump(data, out, cls=moleculeEncoder, indent=2)


def getFirstSynonym(mol_pcp):
    """Returns first synonym available.

    :param mol_pcp: (pubchempy.Compound) Molecule from PubChem.
    :return:
        name: (str) Synonym.
    """

    name = mol_pcp.synonyms
    if len(name) != 0:
        name = name[0]
        name = name.lower()
    else:
        name = 'TODO'
    return name


def compare(mol, mol_pcp):
    """Compares formulas of two molecules.

    :param mol: (molecule.Molecule) Molecule from fieFile.
    :param mol_pcp: (pubchempy.Compound) Molecule from PubChem.
    :return:
    """

    cid_pcp = mol_pcp.cid
    inchi_pcp = mol_pcp.inchikey

    name_pcp = mol_pcp.iupac_name
    if not name_pcp:
        name_pcp = getFirstSynonym(mol_pcp)

    form_pcp = mol_pcp.molecular_formula

    # print(mol_pcp.cid)
    smiles_match = match_smiles(mol.smiles, mol_pcp.canonical_smiles)
    if not smiles_match:
        return False

    formula_match = match_formula(mol.form, form_pcp)
    if not formula_match:
        return False

    mol.cid_pcp = cid_pcp
    mol.inchi_pcp = inchi_pcp
    mol.name_pcp = name_pcp
    mol.form_pcp = form_pcp

    cas = get_cas_pcp(cid_pcp)
    mol.cas = cas

    print('\t{} {} {}'.format(mol.name_pcp, mol.cas, mol_pcp.cid))
    return True


def match_smiles(reference, smiles):
    """Checks if two SMILES strings match.

    :param reference: (str) Reference SMILES string.
    :param smiles: (str) Other SMILES string.
    :return:
    """

    if reference != smiles:
        smiles = utils.getCanonicalSmiles(smiles)
        if reference != smiles:
            print('Smiles do not match: {} != {}'.format(reference, smiles))
            # sys.exit(123)
            return False

    return True


def match_formula(reference, formula):
    """Checks if two molecular formulas match.

    :param reference: (str) Reference formula.
    :param formula: (str) Other formula.
    :return: (boolean)
    """

    s1 = re.findall(r'([A-Z][a-z]*)(\d*)', reference)

    # for alcohol with formula Cx{OH}yHz
    # tmp = re.findall(r'([A-Z][a-z]*)(\d*)', reference)
    # nO = 3
    # s = 'C{}H{}O{}'.format(tmp[0][1], int(tmp[3][1]) + nO * int(tmp[2][1]), nO)
    # s1 = re.findall(r'([A-Z][a-z]*)(\d*)', s)

    s2 = re.findall(r'([A-Z][a-z]*)(\d*)', formula)

    for i in range(len(s1)):
        s1[i] = (s1[i][0].encode(), s1[i][1].encode())
        if not s1[i][1]:
            s1[i] = (s1[i][0], '1'.encode())

    for i in range(len(s2)):
        s2[i] = (s2[i][0].encode(), s2[i][1].encode())
        if not s2[i][1]:
            s2[i] = (s2[i][0], '1'.encode())

    s1 = sorted(s1)
    s2 = sorted(s2)

    if len(s1) != len(s2):
        print('\tFormulas do not match: {} {}'.format(reference, formula))
        return False

    for i in range(len(s1)):
        if s1[i] != s2[i]:
            print(s1[i], s2[i])
            print('\tFormulas do not match: {} {}'.format(reference, formula))
            return False

    return True


def get_cas_pcp(cid):
    """Extracts CAS fom synonyms of pubchempy.Compound or empty string.

    :param cid: (int) CID from PubChem.
    :return:
        cas: (str) CAS Registry Number.
    """

    results = pcp.get_synonyms(cid, 'cid')
    for result in results:
        for syn in result.get('Synonym', []):
            match = re.match('(\d{2,7}-\d\d-\d)', syn)
            if match:
                cas = match.group(1)
                return cas
    return ''


def add_other_names(mol, mol_pcp):
    """Adds other names available in PubChem to molecule.

    :param mol: (molecule.Molecule) Molecule from fieFile.
    :param mol_pcp: (pubchempy.Compound) Molecule from PubChem.
    """

    names = ['allowed_name', 'cas_like_style_name', 'systematic_name', 'traditional_name']
    for n in names:
        try:
            method_to_call = getattr(mol_pcp, n)
            result = method_to_call
            mol.add_other_name(result)
        except AttributeError:
            continue
