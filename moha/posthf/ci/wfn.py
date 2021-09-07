from moha.system.basis import SlaterDeterminant
import numpy as np 

class CIWaveFunction(object):
    """
    """
    def __init__(self,nelec,norb):
        self.nelec = nelec
        self.norb = norb


class CISDWaveFunction(CIWaveFunction):
    """
    """
    def __init__(self,coefficient,ci):
        pass
