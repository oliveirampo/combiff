from abc import ABC
import re


class Family(ABC):
    """Base class for Family object.

    Attributes:
        letter: (str) Letter associated to family.
        cod: (str) Code of family.
        nam: (str) Name of family.
    """

    def __init__(self, letter, cod, nam):
        """Constructs all the necessary attributes for this object.

        :param letter: (str) Letter associated to family.
        :param cod: (str) Code of family.
        :param nam: (str) Name of family.
        """

        self.letter = letter
        self.cod = cod
        self.nam = nam

    def __str__(self):
        return '\n\t{}: {} - {}\n'.format(self.letter, self.cod, self.nam)

    @staticmethod
    def get_num_of_carbons(frm):
        """Returns number of carbon atoms.

        :param frm: (str) Molecular formula.
        :return:
            (int) Number of carbon atoms.
        """

        form = frm.replace('_','')
        atoms = re.findall(r'([A-Z][a-z]?)(\d*)', form)

        nC = 0
        for i in range(len(atoms)):
            atm = atoms[i]

            if atm[0] == 'C':
                nC = atm[1]
                if nC == '':
                    nC = 1
                else:
                    nC = int(nC)

                break

        return nC

    @staticmethod
    def get_num_of_other_atoms(form):
        """Returns number of atoms excluding Carbon atoms.

        :param form: (str) Molecular formula.
        :return:
            (int) Number of atoms.
        """

        atoms = re.findall(r'([A-Z][a-z]?)(\d*)', form)

        # remove Carbon atom
        ignore = []
        for i in range(len(atoms)):
            atm = atoms[i]

            if atm[0] == 'C':
                ignore.append(i)

        count = Family.count_atoms(atoms, ignore)
        return count

    @staticmethod
    def get_num_of_other_atoms_without_H(form):
        """Returns number of atoms excluding Carbon and Hydrogen atoms.

        :param form: (str) Molecular formula.
        :return:
            (int) Number of atoms.
        """

        atoms = re.findall(r'([A-Z][a-z]?)(\d*)', form)

        # remove Carbon and Hydrogen atoms
        ignore = []
        for i in range(len(atoms)):
            atm = atoms[i]

            if atm[0] == 'C':
                ignore.append(i)

            elif atm[0] == 'H':
                ignore.append(i)

        count = Family.count_atoms(atoms, ignore)

        return count

    @staticmethod
    def count_atoms(atoms, ignore):
        count = 0

        for i in range(len(atoms)):
            if i in ignore:
                continue

            elif atoms[i][1] == '':
                count += 1

            else:
                n = atoms[i][1]
                n = int(n)
                count += n

        return count


class Dummy(Family):
    """Implemtation of Dummy Family."""

    def __init__(self):
        """Constructs all the necessary attributes for this object."""

        letter = 'D'
        cod = 'DUMMY'
        nam = 'dummy'

        super(Dummy, self).__init__(letter, cod, nam)


class ALK(Family):
    def __init__(self):
        """Constructs all the necessary attributes for this object."""

        letter = 'A'
        cod = 'ALK'
        nam = 'alkane'

        super(ALK, self).__init__(letter, cod, nam)


class ROH(Family):
    def __init__(self):
        """Constructs all the necessary attributes for this object."""

        letter = 'L'
        cod = 'ROH'
        nam = 'alcohol'

        super(ROH, self).__init__(letter, cod, nam)


class HAL(Family):
    def __init__(self):
        """Constructs all the necessary attributes for this object."""

        letter = 'X'
        cod = 'HAL'
        nam = 'halomethane'

        super(HAL, self).__init__(letter, cod, nam)
