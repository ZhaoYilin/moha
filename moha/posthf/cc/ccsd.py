from moha.posthf.cc.auxiliary import *
from moha.log import log, timer

import numpy as np
import copy

__all__ = ['CCSDSolver']

class CCSDSolver(object):
    """
    """
    def __init__(self,ham,hfwavefunction,maxiter=100,cutoff=1e-6):
        """
        """
        self.ham = ham
        self.hfwavefunction = hfwavefunction
        self.maxiter = maxiter
        self.cutoff = cutoff

    def kernel(self):
        """
        """
        log.hline()
        log('CCSD Calculation'.format())
        log.hline()
        ham = copy.deepcopy(self.ham)
        hfwavefunction = copy.deepcopy(self.hfwavefunction)
        occ = hfwavefunction.occ
        Nelec = occ['alpha'] + occ['beta']
        C = hfwavefunction.coefficient
        dim = ham.dim*2
        #Transfer Fock integral from spatial to spin basis
        fs = spinfock(hfwavefunction.eorbitals)
        #Transfer electron repulsion integral from atomic basis
        #to molecular basis
        ham.operators['electron_repulsion'].basis_transformation(C)
        #build double bar integral <ij||kl>
        spinints = ham.operators['electron_repulsion'].double_bar 
        ts,td,Dai,Dabij = initialize(fs,spinints,Nelec,dim)
        ECCSD = 0
        DECC = 1.0
        Iter = 0
        log.hline()
        log('{0:2s} {1:3s} {2:4s}'.format('Iter', 'ECC', 'Delta'))
        log.hline()
        while DECC > self.cutoff or Iter> self.maxiter: # arbitrary convergence criteria
            Iter += 1
            OLDCC = ECCSD
            Fae,Fmi,Fme,Wmnij,Wabef,Wmbej = updateintermediates(True,Nelec,dim,fs,spinints,ts,td)
            ts = makeT1(True,Nelec,dim,fs,spinints,ts,td,Dai,Fae,Fmi,Fme)
            td = makeT2(True,Nelec,dim,fs,spinints,ts,td,Dabij,Fae,Fmi,Fme,Wmnij,Wabef,Wmbej)
            ECCSD = ccsdenergy(Nelec,dim,fs,spinints,ts,td)
            DECC = abs(ECCSD - OLDCC)
            log('{0:2d} {1:3f} {2:4f}'.format(Iter, ECCSD, DECC))
        log('CC iterations converged'.format())
        log('CCSD energy = {}'.format(ECCSD))
