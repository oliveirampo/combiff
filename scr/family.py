from abc import ABC, abstractmethod
import sys
import re

class Family(ABC):
    def __init__(self):
        self.letter = ''
        self.cod = ''
        self.nam = ''
        self.fdf = ''

    def __str__(self):
        return '\n\t{}: {} - {}\n'.format(self.letter, self.cod, self.nam)

    def get_num_of_carbons(self, frm):
        form = frm.replace('_','')
        atoms = re.findall(r'(\D*)(\d*)', form)

        if atoms[0][0] != 'C':
            if form[1] == 'l':
                sys.exit('ERROR: TODO')

            elif form[1].isalpha():
                return 1

            sys.exit('ERROR: first atom is not C')

        #print(atoms[0][1], frm)
        return int(atoms[0][1])

    @abstractmethod
    def get_num_of_other_atoms(self, form):
        pass

    def get_num_of_other_atoms_without_H(self, form):
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
            #print(atoms[i])
            n = atoms[i][1]
            n = int(n)
            count += n
        #sys.exit('STOP')
        return count


class Dummy(Family):
    def __init__(self):
        self.letter = 'D'
        self.cod = 'DUMMY'
        self.nam = 'dummy'

    def get_num_of_other_atoms(self, form):
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
        self.letter = 'A'
        self.cod = 'ALK'
        self.nam = 'alkanes'

    def get_num_of_other_atoms(self, form):
        return 0


