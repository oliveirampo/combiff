"""Module with exceptions.

Classes:
    ArgError
    NoFile
    WrongNumberofColumns
    EquationNotImplemented
    MethodNotImplemented
    PropertyNotImplemented
    NoKey
    WrongProperty
    VariableNotDefined
"""

import configparser


class ArgError(Exception):
    """Problem with number of arguments in command line."""

    def __init__(self, nArgs, n):
        s = '\n\tWrong number of arguments: it should be {} instead of {}\n\tUsage:'.format(nArgs, n)
        super().__init__(s)


class NoFile(Exception):
    """File not found."""

    def __init__(self, fileName):
        s = '\n\tNo such file:\n\t{}'.format(fileName)
        super().__init__(s)


class WrongNumberofColumns(Exception):
    """Problem with number of columns in file."""

    def __init__(self, fileName):
        s = 'Row with wrong number of columns in {}'.format(fileName)
        super(WrongNumberofColumns, self).__init__(s)


class EquationNotImplemented(Exception):
    """Equation method not implemented."""

    def __init__(self, typ):
        s = 'Equation not implemented: {}'.format(typ)
        super(EquationNotImplemented, self).__init__(s)


class MethodNotImplemented(Exception):
    """Method not implemented."""

    def __init__(self, args):
        args = '\n{}\nImplement the above method in the above class.\n'.format(args)
        super(MethodNotImplemented, self).__init__(args)


class PropertyNotImplemented(Exception):
    """Property not implemented."""

    def __init__(self, args):
        args = '\n{}\nProperty not implemented.\n'.format(args)
        super(PropertyNotImplemented, self).__init__(args)


class NoKey(KeyError):
    """Mapping key not found."""

    def __init__(self, key, src):
        s = 'KeyError: {} not in {}'.format(key, src)
        super(NoKey, self).__init__(s)


class WrongProperty(KeyError):
    """Mapping property not found."""

    def __init__(self, prop, key, function):
        s = 'Could not convert {} to float: {} in {}'.format(prop, key, function)
        super(WrongProperty, self).__init__(s)


class VariableNotDefined(configparser.NoOptionError):
    """Variable not defined in configuration input file."""

    def __init__(self, option):
        section = 'globalVariables'
        super(VariableNotDefined, self).__init__(option, section)
        print('Configuration file has to be fixed.')
