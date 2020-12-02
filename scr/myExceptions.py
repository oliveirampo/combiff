class ArgError(Exception):
    def __init__(self, nArgs, n):
        s = '\n\tWrong number of arguments: it should be {} instead of {}\n\tUsage:'.format(nArgs, n)
        super().__init__(s)

class NoFile(Exception):
    def __init__(self, fileName):
        s = '\n\tNo such file: {}'.format(fileName)
        super().__init__(s)
