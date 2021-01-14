import configparser


class ArgError(Exception):
    def __init__(self, nArgs, n):
        s = '\n\tWrong number of arguments: it should be {} instead of {}\n\tUsage:'.format(nArgs, n)
        super().__init__(s)

class NoFile(Exception):
    def __init__(self, fileName):
        s = '\n\tNo such file:\n\t{}'.format(fileName)
        super().__init__(s)

class WrongNumberofColumns(Exception):
    def __init__(self, fileName):
        s = 'Row with wrong number of columns in {}'.format(fileName)
        super(WrongNumberofColumns, self).__init__(s)

class EquationNotImplemented(Exception):
    def __init__(self, type):
        s = 'Equation not implemented: {}'.format(type)
        super(EquationNotImplemented, self).__init__(s)

class MethodNotImplemented(Exception):
    def __init__(self, args):
        args = '\n{}\nImplement the above method in the above class.\n'.format(args)
        super(MethodNotImplemented, self).__init__(args)

class PropertyNotImplemented(Exception):
    def __init__(self, args):
        args = '\n{}\nProperty not implemented.\n'.format(args)
        super(PropertyNotImplemented, self).__init__(args)

class NoKey(KeyError):
    def __init__(self, key, src):
        s = 'KeyError: {} not in {}'.format(key, src)
        super(NoKey, self).__init__(s)

class WrongProperty(KeyError):
    def __init__(self, prop, key):
        s = 'Could not convert {} to float: {}'.format(prop, key)
        super(WrongProperty, self).__init__(s)

class VariableNotDefined(configparser.NoOptionError):
    def __init__(self, option):
        section = 'globalVariables'
        super(VariableNotDefined, self).__init__(option, section)
        print('Configuration file has to be fixed.')
