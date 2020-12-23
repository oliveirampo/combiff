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
            getIdentifiers.run()

        elif job == '-writeIdentifier':
            writeIdentifiers.run()

        elif job == '-dbsSearch':
            dbsConfigurationFile = 'inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            dbsSearch.run(dbsConfig)

        elif job == '-selectData':
            dbsConfigurationFile = 'inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            selectData.manualSelection(dbsConfig)

        elif job == '-plotData':
            dbsConfigurationFile = 'inp/dbs.conf'
            dbsConfig = dbsSearch.dbsConfiguration(dbsConfigurationFile)
            plotData.plotAll(dbsConfig)

    except (myExceptions.ArgError) as err:
        print(err)
    except (myExceptions.NoFile, myExceptions.WrongNumberofColumns) as err:
        print(err)
    except (myExceptions.NoKey, myExceptions.WrongProperty) as err:
        print(err)
    except (myExceptions.VariableNotDefined, myExceptions.EquationNotImplemented) as err:
        print(err)
    except (myExceptions.MethodNotImplemented) as err:
        print(err)


if __name__ == "__main__":
    main()
