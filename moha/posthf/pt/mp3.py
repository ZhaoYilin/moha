from moha.posthf.pt.auxiliary import *
from moha.log import log, timer

import numpy as np
import copy

__all__ = ['MP3Solver']

class MP3Solver(object):
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
        log('MP3 Calculation'.format())
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
                            Emp2 += Eri[i,a,j,b]*(2*Eri[i,a,j,b]-Eri[i,b,j,a])\
                            /(eorbitals[i] + eorbitals[j] -eorbitals[a] - eorbitals[b])

        elif occ['alpha'] != occ['beta']:
            for spin in C:
                Eri = ham.operators['electron_repulsion'].basis_transformation(C[spin])
                for i in range(occ[spin]):
                    for j in range(occ[spin]):
                        for a in range(occ[spin],ham.dim):
                            for b in range(occ[spin],ham.dim):
                                Emp2 += Eri[i,a,j,b]*(Eri[i,a,j,b]-0.5*Eri[i,b,j,a])\
                                /(eorbitals[spin][i] + eorbitals[spin][j] -eorbitals[spin][a] - eorbitals[spin][b])
        Emp3 = 0.0
        if occ['alpha'] == occ['beta']:
            Eri = ham.operators['electron_repulsion'].double_bar
            for i in range(occ['alpha']):
                for j in range(occ['alpha']):
                    for k in range(occ['alpha']):
                        for l in range(occ['alpha']):
                            for a in range(occ['alpha'],ham.dim):
                                for b in range(occ['alpha'],ham.dim):
                                    Emp3 += (1/8.0)*Eri[i,j,a,b]*Eri[k,l,i,j]*Eri[a,b,k,l]\
                                    /((eorbitals[i] + eorbitals[j] -eorbitals[a] - eorbitals[b])\
                                    *(eorbitals[k] + eorbitals[l] -eorbitals[a] - eorbitals[b]))
            for i in range(occ['alpha']):
                for j in range(occ['alpha']):
                    for a in range(occ['alpha'],ham.dim):
                        for b in range(occ['alpha'],ham.dim):
                            for c in range(occ['alpha'],ham.dim):
                                for d in range(occ['alpha'],ham.dim):
                                    Emp3 += (1/8.0)*Eri[i,j,a,b]*Eri[a,b,c,d]*Eri[c,d,i,j]\
                                    /((eorbitals[i] + eorbitals[j] -eorbitals[a] - eorbitals[b])\
                                    *(eorbitals[i] + eorbitals[j] -eorbitals[c] - eorbitals[d]))
            for i in range(occ['alpha']):
                for j in range(occ['alpha']):
                    for k in range(occ['alpha']):
                        for a in range(occ['alpha'],ham.dim):
                            for b in range(occ['alpha'],ham.dim):
                                for c in range(occ['alpha'],ham.dim):
                                    Emp3 += Eri[i,j,a,b]*Eri[k,b,c,j]*Eri[a,c,i,k]\
                                    /((eorbitals[i] + eorbitals[j] -eorbitals[a] - eorbitals[b])\
                                    *(eorbitals[i] + eorbitals[k] -eorbitals[c] - eorbitals[c]))

        elif occ['alpha'] != occ['beta']:
            print('not implemented yet')
        log.hline()
        log('{0:2s} {1:3f}'.format('Escf', hfwavefunction.Etot))
        log('{0:2s} {1:3f}'.format('Emp2', Emp2))
        log('{0:2s} {1:3f}'.format('Emp3', Emp3))
        log('{0:2s} {1:3f}'.format('Etot', hfwavefunction.Etot+Emp2+Emp3))
        log.hline()

        return Emp3
