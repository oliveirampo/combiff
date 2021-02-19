from abc import ABC, abstractmethod
import sys
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
        atoms = re.findall(r'(\D*)(\d*)', form)

        if atoms[0][0] != 'C':
            if form[1] == 'l':
                sys.exit('ERROR: TODO')

            elif form[1].isalpha():
                return 1

            sys.exit('ERROR: first atom is not C')

        # print(atoms[0][1], frm)
        return int(atoms[0][1])

    @abstractmethod
    def get_num_of_other_atoms(self, form):
        pass

    @staticmethod
    def get_num_of_other_atoms_without_H(form):
        """Returns number of atoms excluding Carbon and Hydrogen atoms.

        :param form: (str) Molecular formula.
        :return:
            (int) Number of atoms.
        """

        atoms = re.findall(r'(\D*)(\d*)', form)

        # remove first and second items
        if atoms[0][0] != 'C':
            sys.exit('ERROR: first atom is not C')

        if atoms[1][0] != 'H':
            sys.exit('ERROR: second atom is not H')

        atoms = atoms[2:]

        # remove last item is empty
        if not atoms[-1][0]:
            atoms = atoms[:-1]

        count = 0
        for i in range(len(atoms)):
            # print(atoms[i])
            n = atoms[i][1]
            n = int(n)
            count += n
        # sys.exit('STOP')
        return count


class Dummy(Family):
    """Implemtation of Dummy Family."""

    def __init__(self):
        """Constructs all the necessary attributes for this object."""

        letter = 'D'
        cod = 'DUMMY'
        nam = 'dummy'

        super(Family, self).__init__(letter, cod, nam)

    def get_num_of_other_atoms(self, form):
        """Returns number of atoms excluding Carbon atoms.

        :param form: (str) Molecular formula.
        :return:
            n: (int) Number of atoms.
        """

        atoms = re.findall(r'([FlrIONPS])(\d*)', form)
        n = 0
        for a in atoms:
            if a[1] == '':
                n += 1
            else:
                n += int(a[1])
        return n


class ALK(Family):
    def __init__(self):
        """Constructs all the necessary attributes for this object."""

        letter = 'A'
        cod = 'ALK'
        nam = 'alkane'

        super(Family, self).__init__(letter, cod, nam)

    def get_num_of_other_atoms(self, form):
        """Returns number of atoms excluding Carbon and Hydrogen atoms.

        :param form: (str) Molecular formula.
        :return:
            n: (int) Number of atoms.
        """

        return 0


class ROH(Family):
    def __init__(self):
        """Constructs all the necessary attributes for this object."""

        letter = 'L'
        cod = 'ROH'
        nam = 'alcohol'

        super(Family, self).__init__(letter, cod, nam)

    def get_num_of_other_atoms(self, form):
        """Returns number of atoms excluding Carbon atoms.

        :param form: (str) Molecular formula.
        :return:
            n: (int) Number of atoms.
        """

        atoms = re.findall(r'([FlrIONPS])(\d*)', form)
        n = 0
        for a in atoms:
            if a[0] == 'O':
                if a[1] == '':
                    n += 3
                else:
                    n += int(a[1])
        return n


class HAL(Family):
    def __init__(self):
        """Constructs all the necessary attributes for this object."""

        letter = 'X'
        cod = 'HAL'
        nam = 'halomethane'

        super(Family, self).__init__(letter, cod, nam)

    def get_num_of_other_atoms(self, form):
        """Returns number of atoms excluding Carbon atoms.

        :param form: (str) Molecular formula.
        :return:
            n: (int) Number of atoms.
        """

        atoms = re.findall(r'([FlrIONPS])(\d*)', form)
        n = 0
        for a in atoms:
            if a[1] == '':
                n += 1
            else:
                n += int(a[1])
        return n




