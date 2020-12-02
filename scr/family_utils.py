import sys

import family

def initFamilies():
    families = [
        family.Dummy(),
        family.ALK(),
    ]
    return families

def getFamily(familyCod):
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