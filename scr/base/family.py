from rdkit import Chem
from abc import ABC
import numpy as np
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

    def get_letter(self, smiles):
        """Returns letter associated to this Family."""
        return self.letter

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


class MIX(Family):
    def __init__(self):
        """Constructs all the necessary attributes for this object."""

        letter = 'G'
        cod = 'MIX'
        nam = 'several_families'

        super(MIX, self).__init__(letter, cod, nam)

    def get_letter(self, smiles):
        """Return letter depending on type of atoms."""

        min_fg_list = ['F', 'Cl', 'Br', 'I',
                       '[OD2]([#6;!$(C=O)])[#6;!$(C=O)]',  # 4  ROR
                       '[CX3H2](=O)',                      # 5  HCOH
                       '[CX3H1](=O)[#6]',                  # 6  RCOH
                       '[#6][CX3](=O)[#6]',                # 7  RCOR
                       '[$([#6][CX3](=O)[OX2H0][#6])]',    # 8  RCOOR
                       '[CX3H1](=O)[OX2H0][#6]',           # 9  HCOOR
                       '[#6;!$(C=O)][OX2H]',               # 10  ROH
                       '[CX3](=O)[OX2H1]',                 # 11  RCOOH
                       '[NX3;H2,H1,H0;!$(NC=O)]',          # 12 RN
                       '[NX3][CX3](=[OX1])[#6]']           # 13 RCON

        min_fg_mask = np.zeros(shape=(len(min_fg_list)))
        import sys
        # sys.exit('STOP')

        mol = Chem.MolFromSmiles(smiles)

        # check if smiles has no listed functional group
        for i in range(len(min_fg_list)):
            smarts = min_fg_list[i]

            fg1 = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(fg1)

            if matches:
                min_fg_mask[i] += 1

        nFG = np.count_nonzero(min_fg_mask)

        # ALK
        if nFG == 0:
            return 'A'

        elif nFG > 1:
            return 'S'

        elif min_fg_mask[0]:
            return 'F'

        elif min_fg_mask[1]:
            return 'C'

        elif min_fg_mask[2]:
            return 'B'

        elif min_fg_mask[3]:
            return 'I'

        elif min_fg_mask[4]:
            return 'O'

        elif min_fg_mask[5]:
            return 'A'

        elif min_fg_mask[6]:
            return 'A'

        elif min_fg_mask[7]:
            return 'K'

        elif min_fg_mask[8]:
            return 'E'

        elif min_fg_mask[9]:
            return 'E'

        elif min_fg_mask[10]:
            return 'L'

        elif min_fg_mask[11]:
            return 'D'

        elif min_fg_mask[12]:
            return 'N'

        elif min_fg_mask[13]:
            return 'M'

        return self.letter
