from moha.posthf.pt.auxiliary import *
from moha.log import log, timer

import numpy as np
import copy

class MP2Solver(object):
    """
    """
    def __init__(self,ham,hfwavefunction):
        """
        """
        self.ham = ham
        self.hfwavefunction = hfwavefunction

    def kernel(self):
        """
        """
        log.hline()
        log('MP2 Calculation'.format())
        log.hline()
        ham = copy.deepcopy(self.ham)
        hfwavefunction = copy.deepcopy(self.hfwavefunction)
        occ = hfwavefunction.occ
        C = hfwavefunction.coefficient
        eorbitals = hfwavefunction.eorbitals
        Emp2 = 0.0
        if occ['alpha'] == occ['beta']:
            Eri = ham.operators['electron_repulsion'].basis_transformation(C)
            for i in range(occ['alpha']):
                for j in range(occ['alpha']):
                    for a in range(occ['alpha'],ham.dim):
                        for b in range(occ['alpha'],ham.dim):
                            Emp2 += Eri[i,a,j,b]*(2*Eri[i,a,j,b]-Eri[i,b,j,a])/(eorbitals[i] + eorbitals[j] -eorbitals[a] - eorbitals[b])

        elif occ['alpha'] != occ['beta']:
            for spin in C:
                Eri = ham.operators['electron_repulsion'].basis_transformation(C[spin])
                for i in range(occ[spin]):
                    for j in range(occ[spin]):
                        for a in range(occ[spin],ham.dim):
                            for b in range(occ[spin],ham.dim):
                                Emp2 += Eri[i,a,j,b]*(Eri[i,a,j,b]-0.5*Eri[i,b,j,a])/(eorbitals[spin][i] + eorbitals[spin][j] -eorbitals[spin][a] - eorbitals[spin][b])
        log.hline()
        log('{0:2s} {1:3f}'.format('Escf', hfwavefunction.Etot))
        log('{0:2s} {1:3f}'.format('Emp2', Emp2))
        log('{0:2s} {1:3f}'.format('Etot', hfwavefunction.Etot+Emp2))
        log.hline()

        return Emp2
