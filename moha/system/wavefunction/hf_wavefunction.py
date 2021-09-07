import numpy as np
from moha.system.wavefunction.base import BaseWaveFunction

class HFWaveFunction(BaseWaveFunction):
    """
    """
    def __init__(self,dim,occ,coefficient,density,fock,eorbitals,Eelec,Etot):
        self.dim = dim
        self.occ = occ
        self.coefficient = coefficient
        self.density = density
        self.fock = fock
        self.eorbitals = eorbitals
        self.Eelec = Eelec
        self.Etot = Etot

    @property
    def configuration(self):
        c = {}
        for spin in self.occ:
            c[spin] = [1]*self.occ[spin] + [0]*(self.dim - self.occ[spin])
        return c
