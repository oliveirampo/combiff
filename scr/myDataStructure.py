from numpy import ndarray
import numpy as np
import json


from . import selectData


class SelectedData:
    """Object that contains selected data.

    Attributes:
        smiles: () SMILES string.
        code: (str) Molecule code.
        mlp: (str) Melting point.
        mlp_src: (str) Source of melting point.
        blp: (str) Boiling point.
        blp_src: (str) Source of boiling point.
        tem_cri: (str) Critical temperature/
        tem_cri_src: (str) Source of critical temperature.
        eps: (str) Permittivity.
        eps_src: (str) Source of permittivity.
        properties: (dict) Dictionary of properties.

    Functions:
        toDict()
    """

    def __init__(self, smiles):
        """Constructs all the necessary attributes for this object.

        :param smiles: (str) SMILES string.
        """

        self.smiles = smiles
        self.code = ''
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
        """Returns dictionary with selected data.

        :return: (dict)
        """

        data = {'__type__': 'SelectedData', 'smiles': self.smiles, 'code': self.code, 'mlp': self.mlp,
                'mlp_src': self.mlp_src, 'blp': self.blp, 'blp_src': self.blp_src, 'tem_cri': self.tem_cri,
                'tem_cri_src': self.tem_cri_src, 'eps': self.eps, 'eps_src': self.eps_src, 'prop': self.properties}

        return data

    @staticmethod
    def getTransitionPoint(data):
        """Returns first value of transition points provided in table.

        :param data: (pandas DataFrame) Table with transition points.
        :return:
            tem: (float, str) Temperature.
            src: (str) Source.
        """

        values = selectData.getFirstValueOfTransitionPoint(data)
        if len(values) == 2:
            tem = values[0]
            src = values[1]
            return tem, src
        return '', ''

    def addMeltingPoint(self, data):
        """Adds melting point to SelectedData object.

        :param data: (pandas DataFrame)
        """

        if self.mlp != '':
            return

        tem, src = self.getTransitionPoint(data)
        self.mlp = tem
        self.mlp_src = src

    def addBoilingPoint(self, data):
        """Adds boiling point to SelectedData object.

        :param data: (pandas DataFrame)
        """

        if self.blp != '':
            return

        tem, src = self.getTransitionPoint(data)
        self.blp = tem
        self.blp_src = src

    def addCriticalTemperature(self, data):
        """Adds critical temperature to SelectedData object.

        :param data: (pandas DataFrame)
        """

        if self.tem_cri != '':
            return

        tem, src = self.getTransitionPoint(data)
        self.tem_cri = tem
        self.tem_cri_src = src

    def addPermittivity(self, data):
        """Adds permittivity to SelectedData object.

        :param data: (pandas DataFrame)
        """

        if self.eps != '':
            return

        tem, src = self.getTransitionPoint(data)
        self.eps = tem
        self.eps_src = src

    def appendProperty(self, prop, mask, x_values, y_values, pre_values, src_values, fid_values, met_values):
        """Appends property to list of properties.

        :param prop: (str) Property code.
        :param mask: (numpy.ndarray of booleans) Array that indicates which coordinates are highlighted.
        :param x_values: (numpy.ndarray) X values.
        :param y_values: (numpy.ndarray) Y values.
        :param pre_values: (numpy.ndarray) List of pressure values.
        :param src_values: (numpy.ndarray) List of LARS codes.
        :param fid_values: (numpy.ndarray) List of annotations.
        :param met_values: (numpy.ndarray) (numpy.ndarray) List of methods.
        :return:
        """

        pre = pre_values[mask]
        tem = x_values[mask]
        val = y_values[mask]
        src = src_values[mask]
        fid = fid_values[mask]
        met = met_values[mask]

        if val.shape[0] != 0:
            if prop in self.properties:
                self.properties[prop].appendData(pre, tem, val, src, fid, met)
            else:
                self.properties[prop] = Property(prop, pre, tem, val, src, fid, met)

    def addProperty(self, prop, property_arg):
        """Adds property to property list.

        :param prop: (str) Property code.
        :param property_arg: (Property)
        """

        if property_arg.val.shape[0] != 0:
            self.properties[prop] = property_arg

    def hasProperty(self, prop):
        """Checks if property is in list of properties.

        :param prop: (str) Property code.
        :return: (boolean)
        """

        return prop in self.properties

    def hasData(self, prop):
        """Checks if property has values associated to it.

        :param prop: (str) Property code.
        :return: (boolean)
        """

        if prop not in self.properties:
            return False

        property_arg = self.properties[prop]
        val = np.asarray(property_arg.val, dtype=np.float32)
        return len(val) != 0

    def getProperty(self, prop):
        """Returns property.

        :param prop: (str) Property code.
        :return: (Property)
        """

        return self.properties[prop]


class Property:
    """Property object.

    Attributes:
        prop: (str) Property code.
        pre: (float) Pressure.
        tem: (float) Temperature.
        val: (float) Value of property.
        src: (str) LARS code.
        fid: (str) Annotation.
        met: (str) Method.
    Functions:
        appendData(pre, tem, val, src, fid, met)
        toDict()
        toString()
    """

    def __init__(self, prop, pre, tem, val, src, fid, met):
        """Constructs all the necessary attributes for this object.

        :param prop: (str) Property code.
        :param pre: (list) Pressure.
        :param tem: (list) Temperature.
        :param val: (list) Value of property.
        :param src: (list) LARS code.
        :param fid: (list) Annotation.
        :param met: (list) Method.
        """

        self.prop = prop
        self.pre = pre
        self.tem = tem
        self.val = val
        self.src = src
        self.fid = fid
        self.met = met

    def appendData(self, pre, tem, val, src, fid, met):
        """Appends data to Property object.

        :param pre: (list) Pressure.
        :param tem: (list) Temperature.
        :param val: (list) Value associated to property.
        :param src: (list) LARS code.
        :param fid: (list) Annotation
        :param met: (list) Method.
        :return:
        """

        for i in range(len(val)):
            self.pre.append(pre[i])
            self.tem.append(tem[i])
            self.val.append(val[i])
            self.src.append(src[i])
            self.fid.append(fid[i])
            self.met.append(met[i])

    def toDict(self):
        """Converts Property object to dictionary.
        This method is used to make this Class serializable.

        :return:
            data: (dict)
        """

        data = {'__type__': 'Property', 'prop': self.prop, 'pre': self.pre, 'tem': self.tem, 'val': self.val,
                'src': self.src, 'fid': self.fid, 'met': self.met}
        return data

    def toString(self):
        """Converts Property to string."""

        s = '{}\n{}\n{}\n{}\n{}'.format(self.prop, self.pre, self.tem, self.val, self.src)
        return s


class SelectedDataEncoder(json.JSONEncoder):
    """JSON <http://json.org> encoder for SelectedData object.
    Returns a serializable object.

    Functions:
        def default(obj)
    """

    def default(self, obj):
        """Returns a serializable object for ``o``, or calls the base implementation
        (to raise a ``TypeError``)."""

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
    """Method to decode a SelectedData object.

    :param obj:
    :return:
        mol: (SelectedData)
    """

    if '__type__' in obj and obj['__type__'] == 'SelectedData':
        smiles = obj['smiles']
        data = SelectedData(smiles)

        code = ''
        if 'code' in obj:
            code = obj['code']

        mlp = obj['mlp']
        mlp_src = obj['mlp_src']
        blp = obj['blp']
        blp_src = obj['blp_src']
        tem_cri = obj['tem_cri']
        tem_cri_src = obj['tem_cri_src']
        eps = obj['eps']
        eps_src = obj['eps_src']

        data.code = code
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
    """Method to decode a Property object.

    :param obj:
    :return:
        mol: (Property)
    """

    if '__type__' in obj and obj['__type__'] == 'Property':
        prop = obj['prop']
        pre = obj['pre']
        tem = obj['tem']
        val = obj['val']
        src = obj['src']
        fid = obj['fid']
        met = obj['met']

        property_arg = Property(prop, pre, tem, val, src, fid, met)
        return property_arg

    else:
        return obj
