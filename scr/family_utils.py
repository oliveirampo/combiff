"""Module for handling Family object.

Methods:
    initFamilies()
    getFamily(familyCod)
"""

import sys

from . import family


def initFamilies():
    """Returns list of predefined families.
    """
    families = [
        family.Dummy(),
        family.ALK(),
        family.ROH(),
        family.HAL(),
    ]
    return families


def getFamily(familyCod):
    """Returns family given code.

    :param familyCod: (str) Family code.
    :return: (Family object)
    """

    families = initFamilies()
    idx = -1
    for i in range(len(families)):
        #print(families[i].cod)
        if families[i].cod == familyCod:
            idx = i
            break

    if idx == -1:
        s = '\n\tERROR: no family found for familyCod = {}\n'.format(familyCod)
        print(s)
        sys.exit(666)

    return families[idx]
