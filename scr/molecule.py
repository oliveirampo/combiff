"""Package that handles Molecule objects.

Classes:
    Molecule
    moleculeEncoder(json.JSONEncoder)

Methods:
    moleculeDecoder(obj)
"""


import json


class Molecule:
    """Molecule object.

    Attributes:
        smiles: (str) SMILES string.
        name: (str) Molecule name.
        code: (str) Molecule code.
        cas: (str) CAS registry number.
        num_frags: (str, int) Number of fragments.
        frags: (list) Array of fragments.
        num_links: (str, int) Number of linkages.
        links: (list) Array of linkages.

        other_names: (list) List of other names associated to molecule.
        form: (str) Molecule formula.
        match: (boolean)
        inchi: (str) Standard InchiKey.
        cid_pcp: (str) CID identifier from PubChem.
        inchi_pcp: (str) InchiKey from PubChem.
        name_pcp: (str) Molecule mame from PubChem.
        form_pcp: (str) Molecule formula from PubChem.
        nC: (str, int) Number of carbon atoms.
        atoms: (list) List of atoms.
        topo_cod: (str) Topology code used in GROMOS building block.
    """

    def __init__(self, smiles, name, code, cas, num_frags, frags, num_links, links):
        """Constructs all the necessary attributes for this object.

        :param smiles: (str) SMILES string.
        :param name: (str) Molecule name.
        :param code: (str) Molecule code.
        :param cas: (str) CAS registry number.
        :param num_frags: (str, int) Number of fragments.
        :param frags: (list) Array of fragments.
        :param num_links: (str, int) Number of linkages.
        :param links: (list) Array of linkages.
        """

        self.smiles = smiles
        self.name = name
        self.enu_name = name
        self.code = code
        self.enu_code = code
        self.cas = cas
        self.num_frags = num_frags
        self.frags = list(frags)
        self.num_links = num_links
        self.links = list(links)

        self.other_names = []
        self.form = ''
        self.match = False
        self.inchi = ''
        self.cid_pcp = ''
        self.inchi_pcp = ''
        self.name_pcp = ''
        self.form_pcp = ''
        self.nC = ''
        self.atoms = []
        self.topo_cod = ''

    def __str__(self):
        s = '{} {}'.format(self.enu_code, self.smiles)
        return s


class moleculeEncoder(json.JSONEncoder):
    """JSON <http://json.org> encoder for Molecule data structure.
        Returns a serializable object.

        Methods:
            def default(obj)
        """

    def default(self, obj):
        """Returns a serializable object for ``o``, or calls the base implementation
        (to raise a ``TypeError``).

        For example, to support arbitrary iterators, you could
        implement default like this::

            def default(self, o):
                try:
                    iterable = iter(o)
                except TypeError:
                    pass
                else:
                    return list(iterable)
                # Let the base class default method raise the TypeError
                return JSONEncoder.default(self, o)
        """

        if isinstance(obj, Molecule):
            return {
                '__type__': 'Molecule',
                'smiles': obj.smiles,
                'name': obj.name,
                'enu_name': obj.enu_name,
                'code': obj.code,
                'enu_code': obj.enu_code,
                'cas': obj.cas,
                'num_frags': obj.num_frags,
                'frags': obj.frags,
                'num_links': obj.num_links,
                'links': obj.links,
                'other_names': obj.other_names,
                'form': obj.form,
                'match': obj.match,
                'inchi': obj.inchi,
                'cid_pcp': obj.cid_pcp,
                'inchi_pcp': obj.inchi_pcp,
                'name_pcp': obj.name_pcp,
                'form_pcp': obj.form_pcp,
                'nC': obj.nC,
                'atoms': obj.atoms,
                'topo_cod': obj.topo_cod
            }
        return json.JSONEncoder.default(self, obj)


def moleculeDecoder(obj):
    """Method to decode a Molecule object.

    :param obj:
    :return:
        mol: (Molecule)
    """

    if '__type__' in obj and obj['__type__'] == 'Molecule':
        mol = Molecule(obj['smiles'], obj['name'], obj['code'], obj['cas'], obj['num_frags'], obj['frags'], obj['num_links'], obj['links'])
        mol.enu_name = obj['enu_name']
        mol.enu_code = obj['enu_code']
        mol.other_names = obj['other_names']
        mol.form = obj['form']
        mol.match = obj['match']
        mol.inchi = obj['inchi']
        mol.cid_pcp = obj['cid_pcp']
        mol.inchi_pcp = obj['inchi_pcp']
        mol.name_pcp = obj['name_pcp']
        mol.form_pcp = obj['form_pcp']
        mol.nC = obj['nC']
        mol.atoms = obj['atoms']
        mol.topo_cod = obj['topo_cod']
        return mol
    else:
        return obj
