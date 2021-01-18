from numpy import ndarray
import numpy as np
import json
import sys


import selectData


class SelectedData():
    def __init__(self, smiles):
        self.smiles = smiles
        self.mlp = ''
        self.mlp_src = ''
        self.blp = ''
        self.blp_src = ''
        self.tem_cri = ''
        self.tem_cri_src = ''
        self.eps = ''
        self.eps_src = ''
        self.properties = {}


    def __str__(self):
        s = '{} {} {} {} {} {} {}'\
            .format(self.smiles, self.mlp, self.mlp_src, self.blp, self.blp_src, self.tem_cri, self.tem_cri_src)

        s = '{}\n{}'.format(s, len(self.properties))

        for key in self.properties:
            prop = self.properties[key]
            propString = prop.toString()
            s = '{}\n{}'.format(s, propString)

        return s


    def toDict(self):
        data = {}
        data ['__type__'] = 'SelectedData'
        data['smiles'] = self.smiles
        data['mlp'] = self.mlp
        data['mlp_src'] = self.mlp_src
        data['blp'] = self.blp
        data['blp_src'] = self.blp_src
        data['tem_cri'] = self.tem_cri
        data['tem_cri_src'] = self.tem_cri_src
        data['eps'] = self.eps
        data['eps_src'] = self.eps_src
        data['prop'] = self.properties

        return data


    def getTransitionPoint(self, data):
        values = selectData.getFirstValueOfTransitionPoint(data)
        if len(values) == 2:
            tem = values[0]
            src = values[1]
            return tem, src
        return '', ''


    def addMeltingPoint(self, data):
        tem, src = self.getTransitionPoint(data)
        self.mlp = tem
        self.mlp_src = src


    def addBoilingPoint(self, data):
        tem, src = self.getTransitionPoint(data)
        self.blp = tem
        self.blp_src = src


    def addCriticalTemperature(self, data):
        tem, src = self.getTransitionPoint(data)
        self.tem_cri = tem
        self.tem_cri_src = src


    def addPermittivity(self, data):
        tem, src = self.getTransitionPoint(data)
        self.eps = tem
        self.eps_src = src


    def appendProperty(self, prop, mask, x_values, y_values, pre_values, src_values, fid_values, met_values):
        pre = pre_values[mask]
        tem = x_values[mask]
        val = y_values[mask]
        src = src_values[mask]
        fid = fid_values[mask]
        met = met_values[mask]

        if prop in self.properties:
            self.properties[prop].appendData(pre, tem, val, src, fid, met)
        else:
            self.properties[prop] = Property(prop, pre, tem, val, src, fid, met)


    def addProperty(self, prop, property):
        self.properties[prop] = property


    def hasProperty(self, prop):
        return prop in self.properties


    def hasData(self, prop):
        if prop not in self.properties:
            return False

        property = self.properties[prop]
        val = np.asarray(property.val, dtype=np.float32)
        return len(val) != 0


    def getProperty(self, prop):
        return self.properties[prop]


class Property():
    def __init__(self, prop, pre, tem, val, src, fid, met):
        self.prop = prop
        self.pre = pre
        self.tem = tem
        self.val = val
        self.src = src
        self.fid = fid
        self.met = met


    def appendData(self, pre, tem, val, src, fid, met):
        for i in range(len(val)):
            self.pre.append(pre[i])
            self.tem.append(tem[i])
            self.val.append(val[i])
            self.src.append(src[i])
            self.fid.append(fid[i])
            self.met.append(met[i])

    def toDict(self):
        data = {}
        data['__type__'] = 'Property'
        data['prop'] = self.prop
        data['pre'] = self.pre
        data['tem'] = self.tem
        data['val'] = self.val
        data['src'] = self.src
        data['fid'] = self.fid
        data['met'] = self.met
        return data

    def toString(self):
        s = '{}\n{}\n{}\n{}\n{}'.format(self.prop, self.pre, self.tem, self.val, self.src)
        return s


class SelectedDataEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, SelectedData):

            dictionary = obj.toDict()
            return dictionary

        elif isinstance(obj, Property):
            dictionary = obj.toDict()
            return dictionary

        elif isinstance(obj, ndarray):
            return obj.tolist()

        else:
            return json.JSONEncoder.default(self, obj)


def selectedDataDecoder(obj):
    if '__type__' in obj and obj['__type__'] == 'SelectedData':
        smiles = obj['smiles']
        data = SelectedData(smiles)

        mlp = obj['mlp']
        mlp_src = obj['mlp_src']
        blp = obj['blp']
        blp_src = obj['blp_src']
        tem_cri = obj['tem_cri']
        tem_cri_src = obj['tem_cri_src']
        eps = obj['eps']
        eps_src = obj['eps_src']

        data.mlp = mlp
        data.mlp_src = mlp_src
        data.blp = blp
        data.blp_src = blp_src
        data.tem_cri = tem_cri
        data.tem_cri_src = tem_cri_src
        data.eps = eps
        data.eps_src = eps_src

        propertiesDict = obj.get('prop')
        propertiesString = json.dumps(propertiesDict)

        propertiesDict = json.loads(propertiesString, object_hook=propertyDecoder)
        data.properties = propertiesDict

        return data

    else:
        return obj


def propertyDecoder(obj):
    if '__type__' in obj and obj['__type__'] == 'Property':
        prop = obj['prop']
        pre = obj['pre']
        tem = obj['tem']
        val = obj['val']
        src = obj['src']
        fid = obj['fid']
        met = obj['met']

        property = Property(prop, pre, tem, val, src, fid, met)
        return property

    else:
        return obj
