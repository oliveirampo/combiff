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

from . import myExceptions
from . import getIdentifiers
from . import writeIdentifiers
from . import dbsSearch
from . import selectData
from . import plotData
from . import IO


def main():
    """Main function of pipeline to extract data for list of molecules using SMILES strings."""

    try:
        if len(sys.argv) == 1:
            raise myExceptions.ArgError(2, len(sys.argv))

        job = sys.argv[1]

        if job == '-getIdentifier':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            getIdentifiers.run(dbsConfig)

        elif job == '-writeIdentifier':
            writeIdentifiers.run()

        elif job == '-searchProp':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            dbsSearch.run(dbsConfig)

        elif job == '-selectData':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            selectData.manualSelection(dbsConfig)

        elif job == '-addSelectedData':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            selectData.getSelectedData(dbsConfig)

        elif job == '-addSelectedDataForProperty':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            selectData.getSelectedDataForProperty(dbsConfig)

        elif job == '-assignCode':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            selectData.assignCode(dbsConfig)

        elif job == '-plotData':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            plotData.plotAll(dbsConfig)

        elif job == '-writeOutputFiles':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)

            IO.writeMolDataFile(dbsConfig)

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
