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

        elif job == '-plotData':
            dbsConfigurationFile = '../inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            plotData.plotAll(dbsConfig)

    # except (myExceptions.MethodNotImplemented) as err:
    #     print(err)
    except (myExceptions.VariableNotDefined) as err:
        print(err)
    except (KeyError, Exception) as err:
        print(err)


if __name__ == "__main__":
    main()
