"""The IO module provides methods for general input and output handling.

Methods:
    readFieFile(fileName)
    readFlsFile(fileName, isomers)
    readFlsFile_helper(pos, data)
    readCiDSmilesFile(fileName)
    writeSelectedDataToJson(dbsConfig, selectedData)
    openSelectedDataFile(dbsConfig)
    readCodSmilesMap(fileName)
    writeMolDataFile(dbsConfig)
"""

import pandas as pd
import numpy as np
import json
import math
import sys
import os

from scr.base.molecule import Molecule
from scr.base import myDataStructure
from scr.base import myExceptions
from scr.base import utils
import xml.etree.ElementTree as ET


def readXmlFile(fileName):
    """Reads family isomer enumeration (fie) file.

    :param fileName: (str) File name.
    :return:
        df: (pandas DataFrame) Table with ENU code, molecular formula and SMILEs string.
    """

    tree = ET.parse(fileName)
    root = tree.getroot()

    data = [[isomer.get('isomer_id'), isomer_list.get('formula'), isomer.find('constitutional_SMILES').text]
            for isomer_lists in root.findall('isomer_lists')
            for isomer_list in isomer_lists.findall('isomer_list')
            for isomer in isomer_list.findall('isomer')]

    df = pd.DataFrame(data, columns=['nam', 'frm', 'smiles'])

    return df


def readFieFile(fileName):
    """Reads family isomer enumeration (fie) file.

    :param fileName: (str) File name.
    :return:
        df: (pandas DataFrame) Table with ENU code, molecular formula and SMILEs string.
    """

    if not os.path.exists(fileName):
        raise myExceptions.NoFile(fileName)
    # df = pd.read_csv(fileName, sep='\s+', comment='#', names=['nam', 'frm', 'smiles'])
    df = pd.read_csv(fileName, sep='\s+', names=['nam', 'frm', 'smiles'])

    df = df[~df['nam'].str.startswith('#')]
    df.reset_index(inplace=True, drop=True)

    return df


def readFlsFile_old(fileName, isomers):
    """Reads molecule file (fls file).

    :param fileName: (str) File name.
    :param isomers: (pandas Dataframe) Table with ENU code, molecular formula and SMILEs string..
    :return:
        molecules: (list) List of molecules.
    """

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


def readFlsFile(fileName):
    """Reads molecule file (fls file).

    :param fileName: (str) File name.
    :return:
        molecules: (dict) Map of smiles to molecule.
    """

    with open(fileName, 'r') as f:
        lines = f.readlines()

    lines = [row.strip().split() for row in lines]

    molecules = {}
    for i in range(len(lines)):
        if not lines[i]:
            continue

        if lines[i][0] == '#':

            # couldn't find a decomposition
            if len(lines[i]) > 2 and lines[i][1] == "couldn\'t":
                smiles = lines[i][6]
                mol = Molecule(smiles, '', '', '', 0, [], 0, '')
                molecules[smiles] = mol

            else:
                continue

        if lines[i][0][0] == '#':
            continue

        if lines[i][0] == 'NEWMOLECULE':
            pos = i + 1
            mol = readFlsFile_helper(pos, lines)

            molecules[mol.smiles] = mol

    return molecules


def readFlsFile_helper(pos, data):
    """Helper function to read fls file.

    :param pos: (int) Current position in the while-loop.
    :param data: (list) Rows with data.
    :return:
    """

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
            frags.append([idx, frg])
        elif data[i][0] == 'NUMLINKS':
            num_links = data[i][1]
        elif data[i][0] == 'LINK':
            frg_1 = data[i][1]
            frg_2 = data[i][2]
            lnk = data[i][3]
            links.append([frg_1, frg_2, lnk])


def read_mtb_file(fileName):
    """Reads MTB file.

    :param fileName: (str) File name.
    :return:
        molecules: (dict) Map of smiles to molecule
    """

    if not os.path.exists(fileName):
        sys.exit(f'No such file: {fileName}')

    with open(fileName, 'r') as f:
        lines = f.readlines()

    i = 0
    molecules = {}
    while i < len(lines):
        row = lines[i]

        if row.startswith('TITLE'):
            title_block = lines[:9]
            i = 9
            if title_block[-1].strip() != 'END':
                sys.exit('TITLE block is not correct.')

        elif row.startswith('MTBUILDBLSOLUTE'):
            i, smiles, mol_block = read_mtb_file_helper(i, lines)
            molecules[smiles] = mol_block

        else:
            sys.exit('Beginning of row not supported: {}'.format(row))

    return title_block, molecules


def read_mtb_file_helper(i, lines):
    j = i
    data = []

    smiles = lines[i + 1].split()[-1]
    if not lines[i + 1].startswith('# smiles'):
        smiles = lines[i + 2].strip()

    while j < len(lines):
        row = lines[j]

        if row.startswith('END'):
            i = j + 1
            j = len(lines)

        data.append(row)
        j += 1

    return i, smiles, data


def readCiDSmilesFile(fileName):
    """Reads file that maps PubChem identifiers (CID) to SMILES string.

    :param fileName: (str) File name.
    :return:
        (pandas DataFrame) Table with CIDs and SMILES strings.
    """

    print('Reading PubChem file...')
    df = pd.read_csv(fileName, sep='\s+', header=None, names=['cid', 'smiles'], dtype={'cid': 'Int64', 'smiles': 'str'})

    if 'canonic' in fileName:
        return df

    print('Canonicalizing SMILES')
    df['smiles'] = df['smiles'].apply(utils.getCanonicalSmiles)

    return df


def writeSelectedDataToJson(dbsConfig, selectedData):
    """Writes selected data to json file.

    :param dbsConfig: (dbsConfiguration object) DBS configuration object.
    :param selectedData: (dict of selectedData) Dictionary that maps SMILES to object that contains the selected data.
    """

    fileName = dbsConfig.getOutFileName('molJsonFile')

    with open(fileName, "w") as outfile:
        json.dump(selectedData, outfile, cls=myDataStructure.SelectedDataEncoder, indent=1)


def openSelectedDataFile(dbsConfig):
    """Loads json file to selectedData object.

    :param dbsConfig:
    :return:
        selectedData: (dict of selectedData) Dictionary that maps SMILES to object that contains the selected data
    """

    fileName = dbsConfig.getOutFileName('molJsonFile')
    if not os.path.exists(fileName):
        return {}

    with open(fileName) as jsonFile:
        selectedData = json.load(jsonFile, object_hook=myDataStructure.selectedDataDecoder)

    return selectedData


def readCidSmilesMap(fileName, columns):
    """Reads plain file that maps SMILES to molecule code.
    # Checks if column of SMILES is unique.
    # Checks if column of codes is unique.

    :param fileName: (str) Name of file.
    :param columns: (arr) Name of columns.
    :return:
    """

    df = pd.read_csv(fileName, sep='\s+', names=columns)

    # check if columns only have unique values.
    isUnique = df['code'].is_unique
    if isUnique is False:
        print('\n\tCode column has duplicate values in {}'.format(fileName))
        sys.exit(123)
    isUnique = df['smiles'].is_unique
    if isUnique is False:
        print('\n\tSMILES column has duplicate values in {}'.format(fileName))
        sys.exit(123)

    return df


def writeMolDataFile(dbsConfig):
    """Writes selected data from json file to plain text file.
    The files created are:
        - mol file to run SAMOS or GROMOS simulation package.
        - file with source code for each property.

    :param dbsConfig: (dbsConfiguration object) DBS configuration object.

    prop == 'dns' gives the same results as prop == 'hvp'
    """

    allSelectedData = openSelectedDataFile(dbsConfig)

    propertyList = dbsConfig.getPropListToBeWrittenToFile()

    dnsVariable = dbsConfig.getVariable('density')
    hvpVariable = dbsConfig.getVariable('vaporizationEnthalpy')

    # write dns, hvp, mlp, blp, tem_cri and eps for SAMOS simulation
    writeMolDat = False
    # print(propertyList)
    for prop in propertyList:
        if (prop == dnsVariable) or (prop == hvpVariable):
            writeMolDat = True
            break

    # [code, smiles, run, pre, tem, wDns, dns, wHvp, hvp, mlp, blp, tem_cri, eps]
    outputData = []
    if writeMolDat:
        writeDnsHvpData(dbsConfig, allSelectedData)

    # write dns, hvp, eps, and prop for GROMOS simulation.
    for prop in propertyList:
        expFileName = 'out/mol_{}.exp'.format(prop)
        srcFileName = 'out/mol_{}.src'.format(prop)
        with open(expFileName, 'w') as expOut, open(srcFileName, 'w') as srcOut:

            for smiles in allSelectedData:
                selectedData = allSelectedData[smiles]
                code = selectedData.code

                properties = selectedData.properties

                dnsVal = 0.0
                if dnsVariable in properties:
                    dnsVal = properties[dnsVariable].val[0]

                nC = int(code[1]) + int(code[2])
                kappa = '4.57E-04'
                eps = selectedData.eps
                epsSrc = selectedData.eps_src
                eps, epsSrc = toString(eps, epsSrc, 'TODO')

                if prop not in properties:
                    continue

                propPre = properties[prop].pre[0]
                propTem = properties[prop].tem[0]
                propVal = properties[prop].val[0]
                propSrc = properties[prop].src[0]

                dnsVal = float(dnsVal)

                # out.write('{:6} {:2} {:6} {:6} {:10} {:8} {:8.2f} {:8} {}\n'
                #           .format(code, nC, propPre, propTem, eps, kappa, dnsVal, propVal, smiles))
                expOut.write('{:6} {:8.2f} {:2} {:10} {:10} {:6} {:6} {:8} {}\n'
                             .format(code, dnsVal, nC, eps, kappa, propTem, propPre, propVal, smiles))

                srcOut.write('{:6} {:6} {:6} {:8} {}\n'
                             .format(code, propTem, propPre, propVal, propSrc))


def organizeValues(prop1, prop2):
    """Makes pairs of values if they have close values of pressure and temperature.
    Pressure tolerance = 5.0 kPa.
    Temperature tolerance = 2.0 K.

    :param prop1: (Property)
    :param prop2: (Property)
    :return:
    """

    preProp1 = prop1.pre
    temProp1 = prop1.tem
    valProp1 = prop1.val

    preProp2 = prop2.pre
    temProp2 = prop2.tem
    valProp2 = prop2.val

    nProp1 = len(valProp1)
    nProp2 = len(valProp2)

    pairIdx = np.full((nProp1, 1), -1)

    i = 0
    while i < nProp1:

        j = 0
        while j < nProp2:

            isPreClose = math.isclose(preProp1[i], preProp2[j], abs_tol=5.0)  # kPa
            isTemClose = math.isclose(temProp1[i], temProp2[j], abs_tol=2.0)  # K
            if isPreClose:

                if isTemClose:
                    if j not in pairIdx[:, 0]:
                        pairIdx[i] = j
                        j = nProp2
                        i += i
            j += 1
        i += 1

    return pairIdx


def writeDnsHvpData(dbsConfig, allSelectedData):
    """Writes density and vaporization enthalpy values to plain file.

    :param dbsConfig: (dbsConfiguration object) DBS configuration object.
    :param allSelectedData: (dict of selectedData) Maps SMILES to object that contains the selected data.
    """

    outFileName = dbsConfig.getOutFileName('molDataFile')
    srcFileName = dbsConfig.getOutFileName('molSrcFile')

    with open(outFileName, 'w') as outFile, open(srcFileName, 'w') as srcFile:
        for smiles in allSelectedData:
            selectedData = allSelectedData[smiles]
            code = selectedData.code
            # print('{:5} {}'.format(code, smiles))

            blp = selectedData.blp
            blpSrc = selectedData.blp_src
            mlp = selectedData.mlp
            mlpSrc = selectedData.mlp_src
            tem_cri = selectedData.tem_cri
            tem_criSrc = selectedData.tem_cri_src
            eps = selectedData.eps
            epsSrc = selectedData.eps_src
            properties = selectedData.properties

            mlp, mlpSrc = toString(mlp, mlpSrc, 0.0)
            blp, blpSrc = toString(blp, blpSrc, 0.0)
            tem_cri, tem_criSrc = toString(tem_cri, tem_criSrc, 0.0)
            eps, epsSrc = toString(eps, epsSrc, 'TODO')

            dnsVariable = dbsConfig.getVariable('density')
            hvpVariable = dbsConfig.getVariable('vaporizationEnthalpy')

            dns = myDataStructure.Property('dns', [], [], [], [], [], [])
            hvp = myDataStructure.Property('hvp', [], [], [], [], [], [])
            if dnsVariable in properties:
                dns = properties[dnsVariable]
            if hvpVariable in properties:
                hvp = properties[hvpVariable]

            jobLetter = ord('a')
            pairIdx = organizeValues(dns, hvp)

            # write dns data
            for i in range(pairIdx.shape[0]):
                hvpIdx = pairIdx[i][0]

                preDns = dns.pre[i]
                temDns = dns.tem[i]
                valDns = dns.val[i]
                srcDns = dns.src[i]

                # write dns/hvp pair if exists
                if hvpIdx != -1:
                    valHvp = hvp.val[hvpIdx]
                    srcHvp = hvp.src[hvpIdx]

                    writeDnsHvpDataHelper(code, jobLetter, smiles, preDns, temDns, 1.0, valDns, srcDns, 1.0,
                                          valHvp, srcHvp, mlp, mlpSrc, blp, blpSrc, tem_cri, tem_criSrc, eps, epsSrc,
                                          outFile, srcFile)

                # write the remaining density values
                else:
                    valHvp = 0.0
                    srcHvp = '%'

                    writeDnsHvpDataHelper(code, jobLetter, smiles, preDns, temDns, 1.0, valDns, srcDns, 0.0,
                                          valHvp, srcHvp, mlp, mlpSrc, blp, blpSrc, tem_cri, tem_criSrc, eps, epsSrc,
                                          outFile, srcFile)

                jobLetter += 1

            # write hvp data:
            for i in range(len(hvp.val)):
                if i in pairIdx[:, 0]:
                    continue

                valHvp = hvp.val[i]
                srcHvp = hvp.src[i]

                preHvp = hvp.pre[i]
                temHvp = hvp.tem[i]
                valDns = 0.0
                srcDns = '%'

                writeDnsHvpDataHelper(code, jobLetter, smiles, preHvp, temHvp, 0.0, valDns, srcDns, 1.0, valHvp, srcHvp,
                                      mlp, mlpSrc, blp, blpSrc, tem_cri, tem_criSrc, eps, epsSrc, outFile, srcFile)

                jobLetter += 1


def toString(val, src, replacement):
    """Converts value to string.

    :param val: (str, float) Value to be converted.
    :param src: (str) Reference source.
    :param replacement: (str) String to be used in case value is missing.
    :return:
    """

    try:
        val = float(val)
        val = '{:8.1f}'.format(val)
    except ValueError:
        val = '{:8}'.format(replacement)
        src = '%'

    return val, src


def writeDnsHvpDataHelper(code, jobLetter, smiles, preDns, temDns, runDns, valDns, srcDns, runHvp, valHvp, srcHvp,
                          mlp, mlpSrc, blp, blpSrc, tem_cri, tem_criSrc, eps, epsSrc, outFile, srcFile):
    """Help function to write dns and hvp reference data.

    :param code: (str) Molecule code.
    :param jobLetter: (str) Letter do differentiate the same molecule in different states.
    :param smiles: (str) SMILES string.
    :param preDns: (str, float) Pressure.
    :param temDns: (str, float) Temperature.
    :param runDns: (str, float) Flag to indicate if density should be included in the optimization.
    :param valDns: (str, float) Density value.
    :param srcDns: (str) Source of density.
    :param runHvp: (str, float) Flag to indicate if vaporization enthalpy should be included in the optimization.
    :param valHvp: (str, float) Vaporization enthalpy value.
    :param srcHvp: (str) Source of vaporization enthalpy.
    :param mlp: (str, float) Melting point.
    :param mlpSrc: (str) Source of melting point.
    :param blp: (str, float) Boiling point.
    :param blpSrc: (str) Source of boiling point.
    :param tem_cri: (str, float) Critical temperature.
    :param tem_criSrc: (str) Source of critical temperature.
    :param eps: (str, float) Permittivity.
    :param epsSrc: (str) Source of permittivity.
    :param outFile: (str) Output file with data.
    :param srcFile: (str) Output file with source.
    """

    mlp = float(mlp)
    blp = float(blp)
    tem_cri = float(tem_cri)
    valDns = float(valDns)
    valHvp = float(valHvp)

    outFile.write('{:5}{} {:18} {:3} {:8.3f} {:6.2f} {:3} {:8.2f} {:3} {:8.2f} {:6} {:6} {} {:>6}\n'
                  .format(code, chr(jobLetter), smiles, 1.0, preDns, temDns, runDns, valDns, runHvp, valHvp,
                          mlp, blp, tem_cri, eps))

    srcFile.write('{:5}{} {:18} {:3} {:8.3f} {:6.2f} {:3} {:8.2f} {:8} {:3} {:8.2f} {:8} {:6} {:8} {:6} '
                  '{:8} {:6} {:8} {} {:8}\n'
                  .format(code, chr(jobLetter), smiles, 1.0, preDns, temDns, runDns, valDns, srcDns, runHvp,
                          valHvp, srcHvp, mlp, mlpSrc, blp, blpSrc, tem_cri, tem_criSrc, eps, epsSrc))


def write_selected_molecule_file(dbsConfig):
    inpFileName = dbsConfig.getOutFileName('inpMoleculeFile')
    outFileName = dbsConfig.getOutFileName('outMoleculeFile')
    fileExtension = outFileName.split('.')[-1]

    if fileExtension == 'fls':
        molecules = readFlsFile(inpFileName)

    elif fileExtension == 'mtb':
        title_block, molecules = read_mtb_file(inpFileName)

    else:
        sys.exit('File extension not supported: {}'.format(fileExtension))

    # get selected molecules
    selected_molecules = {}

    allSelectedData = openSelectedDataFile(dbsConfig)
    for smiles in allSelectedData:
        selected_molecules[smiles] = molecules[smiles]

    if fileExtension == 'fls':
        write_fls(selected_molecules, outFileName)
    elif fileExtension == 'mtb':
        write_mtb_file(title_block, selected_molecules, outFileName)


def write_fls(data, file_name):
    """Write molecule file (file.mtb)"""
    with open(file_name, 'w') as out:
        for smiles in data:
            mol = data[smiles]

            if mol.num_frags == 0:
                out.write("# couldn't find a decomposition for {} with the given fragments\n#\n".format(smiles))

            else:
                out.write('#\nNEWMOLECULE\n#\n')

                out.write('\tNAME\t{}\n'.format(mol.smiles))
                out.write('\tCODE\t{}\n'.format(mol.code))
                out.write('\tCAS\t******\n'.format(mol.cas))

                out.write('\tNUMFRAGS\t{}\n'.format(mol.num_frags))
                out.write('\t#\tidx\tfrg\n')
                for idx, frag in mol.frags:
                    out.write('\tFRAG\t{}\t{}\n'.format(idx, frag))

                out.write('\tNUMLINKS\t{}\n'.format(mol.num_links))
                if mol.num_links != 0:
                    out.write('\t#\tfrg1\tfrg2\tlnkatms\n')
                    for frg_1, frg_2, lnk in mol.links:
                        out.write('\tLINK\t{}\t{}\t{}\n'.format(frg_1, frg_2, lnk))

                out.write('#\nENDMOLECULE\n#\n')


def write_mtb_file(title_block, data, file_name):
    """Write mtb file

    :param title_block: (arr) Title block of MTB file.
    :param data: (arr) MTBUILDBLSOLUTE blocks of data
    :param file_name: (str)
    """

    with open(file_name, 'w') as out:
        for row in title_block:
            out.write(row)

        for smiles in data:
            mol_block = data[smiles]
            for row in mol_block:
                out.write(row)


def update_code_in_mtb_file(dbsConfig):
    allSelectedData = openSelectedDataFile(dbsConfig)
    output_file_name = 'out/mol_CG.mtb'

    # Get building blocks
    title_block, molecules_block = read_mtb_file(output_file_name)
    for smiles in molecules_block:
        selectedData = allSelectedData[smiles]

        block = molecules_block[smiles]
        i = 0
        while i < len(block):
            row = block[i]
            if row.startswith('# RNME'):
                block[i + 1] = selectedData.code + '\n'
                i = len(block)

            i += 1

    write_mtb_file(title_block, molecules_block, output_file_name)
