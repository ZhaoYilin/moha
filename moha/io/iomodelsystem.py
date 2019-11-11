'''
from moha.modelsystem.lattice import Lattice
from moha.modelsystem import sites
class IOModelSystem(object):
    def __init__(self):
        pass

    @classmethod
    def read(cls,latticetype,sitestype,size):
        lattice = getattr(Lattice,latticetype)(size)
        site = getattr(sites,sitestype)

        return lattice,site
'''
