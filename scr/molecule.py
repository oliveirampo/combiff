import json
import sys


class Molecule():
    def __init__(self, smiles, name, code, cas, num_frags, frags, num_links, links):
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
    def default(self, obj):
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
#


