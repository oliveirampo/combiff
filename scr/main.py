import sys

import myExceptions
import getIdentifiers
import writeIdentifiers
import dbsSearch
import selectData
import plotData


def main():
    try:
        if len(sys.argv) == 1:
            raise myExceptions.ArgError(2, len(sys.argv))

        job = sys.argv[1]

        '''
        Get identifiers (name, cas, standard inchikey) from pubchem given smiles strings.
        - If a file with CID and smiles from pubchem is given,
        then the CID is used directly to retrieve the molecule from pubchem.
        - If an empty file is given,
        then the smiles is used to search for the molecule in the database.
        In both cases RDKit is used to canonicalize the smiles strings.
        '''
        if job == '-getIdentifier':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            getIdentifiers.run(dbsConfig)

        # Write identifiers to txt file
        elif job == '-writeIdentifier':
            writeIdentifiers.run()

        '''
        Search for reference data in DBS tables.
        The DBS tables are saved in tmp/ directory.
        The output matched data is saved in data/ directory. 
        '''
        elif job == '-searchProp':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            dbsSearch.run(dbsConfig)

        '''
        Plot all available data for each property and each molecule
        so that the user can select the data points that will be used as reference data.
        '''
        elif job == '-selectData':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            selectData.manualSelection(dbsConfig)

        # Add selected dns and hvp values from files of older version of this script.
        elif job == '-addSelectedData':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            selectData.getSelectedData(dbsConfig)

        # Add additional values of properties from files of older version of this script.
        elif job == '-addSelectedDataForProperty':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            selectData.getSelectedDataForProperty(dbsConfig)

        # Save plot with all data and simulated results if input file is given.
        elif job == '-plotData':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            plotData.plotAll(dbsConfig)

    # except (myExceptions.MethodNotImplemented) as err:
    #     print(err); sys.exit(1)
    except (myExceptions.PropertyNotImplemented) as err:
        print(err); sys.exit(1)
    except (myExceptions.VariableNotDefined) as err:
        print(err); sys.exit(1)
    except (myExceptions.NoKey, myExceptions.WrongProperty) as err:
        print(err); sys.exit(1)


if __name__ == "__main__":
    main()
