"""Pipeline to extract data from different sources. The options are:

-getIdentifier
    Gets identifiers (name, cas, standard inchikey) from pubchem given smiles strings.
        - If a file with CID and smiles from pubchem is given,
        then the CID is used directly to retrieve the molecule from pubchem.
        - If an empty file is given,
        then the smiles is used to search for the molecule in the database.
    In both cases RDKit is used to canonicalize the smiles strings.
-writeIdentifier
    Writes identifiers to txt file
-searchProp
    Searches for reference data in DBS tables.
    The DBS tables are saved in tmp/ directory.
    The output matched data is saved in data/ directory.
-selectData
    Plots all available data for each property and each molecule
    so that the user can select the data points that will be used as reference data.
-addSelectedData
    Adds selected dns and hvp values from files of older version of this script.
-addSelectedDataForProperty
    Adds additional values of properties from files of older version of this script.
-assignCode
    Assigns code to molecule. This code should be easier to grep than the SMILES string.
-plotData
    Saves plot with all data and simulated results if input file is given.
-writeOutputFiles
    Writes selected data from json file to plain text file.

The tables with data from different sources are found in 'wrkDir/tmp/' directory.
"""


import sys

from scr.base import myExceptions
from scr.base import getIdentifiers
from scr.base import writeIdentifiers
from scr.base import dbsSearch
from scr.base import selectData
from scr.base import plotData
from scr.base import IO


def main():
    """Main function of pipeline to extract data for list of molecules using SMILES strings."""

    try:
        if len(sys.argv) == 1:
            raise myExceptions.ArgError(2, len(sys.argv))

        job = sys.argv[1]
        dbsConfigurationFile = '../inp/dbs.conf'
        dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)

        if job == '-getIdentifier':
            getIdentifiers.run(dbsConfig)

        elif job == '-writeIdentifier':
            writeIdentifiers.run()

        elif job == '-searchProp':
            dbsSearch.run(dbsConfig)

        elif job == '-searchPropAlternative':
            alternative_directory = 'tmp2'
            dbsSearch.runAlternative(dbsConfig, alternative_directory)

        elif job == '-selectData':
            selectData.manualSelection(dbsConfig)

        elif job == '-addSelectedData':
            selectData.getSelectedData(dbsConfig)

        elif job == '-addSelectedDataForProperty':
            selectData.getSelectedDataForProperty(dbsConfig)

        elif job == '-assignCode':
            selectData.assignCode(dbsConfig)

        elif job == '-plotData':
            plotData.plotAll(dbsConfig)

        elif job == '-writeOutputFiles':
            IO.writeMolDataFile(dbsConfig)

        elif job == '-saveMoleculeFile':
            IO.write_selected_molecule_file(dbsConfig)

        elif job == '-updateCodeInMTBFile':
            IO.update_code_in_mtb_file(dbsConfig)

        else:
            print('Option not implemented: {}'.format(job))

    except (myExceptions.MethodNotImplemented, myExceptions.PropertyNotImplemented) as err:
        print(err)
        sys.exit(1)
    except myExceptions.NoFile as err:
        print(err)
        sys.exit(1)
    except (myExceptions.VariableNotDefined, myExceptions.NoKey, myExceptions.WrongProperty) as err:
        print(err)
        sys.exit(1)


if __name__ == "__main__":
    main()
