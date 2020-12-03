import sys

import myExceptions
import getIdentifiers
import writeIdentifiers
import dbsSearch
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
            dbsSearch.run()

        elif job == '-plotAll':
            plotData.plotAll()


    except (myExceptions.ArgError) as err:
        print(err)
    except (myExceptions.NoFile) as err:
        print(err)


if __name__ == "__main__":
    main()
