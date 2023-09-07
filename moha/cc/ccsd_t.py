from moha.cc.auxiliary import *
from moha.io.log import log, timer

import numpy as np
import copy

__all__ = []

def ccsd_t(ham,wfn,maxiter=100,E_conv=1.0E-8):
    """Kernel of the solver.

    Returns
    -------
    results : dict
        CCSD calculation results.
    """
    log.hline(char='=')
    log('CCSD(T) Section'.format())
    log.blank()

    nspatial = ham.nspatial
    nspin = nspatial*2
    nelec = wfn.symmetry.n

    C = wfn.coefficient

    #Transfer Fock integral from spatial to spin basis
    fs = spinfock(wfn.orbital_energies)
    #Transfer electron repulsion integral from atomic basis
    #to molecular basis
    ham['Eri'].basis_transformation(C)
    #build double bar integral <ij||kl>
    spinints = ham['Eri'].double_bar 
    ts,td,Dai,Dabij = initialize(fs,spinints,nelec,nspin)
    ECCSD = 0
    DECC = 1.0
    Iter = 0
    log.hline()
    log('{0:2s} {1:3s} {2:4s}'.format('Iter', 'ECC', 'Delta'))
    log.hline()
    while DECC > E_conv or Iter> maxiter: # arbitrary convergence criteria
        Iter += 1
        OLDCC = ECCSD
        Fae,Fmi,Fme,Wmnij,Wabef,Wmbej = updateintermediates(True,nelec,nspin,fs,spinints,ts,td)
        ts = makeT1(True,nelec,nspin,fs,spinints,ts,td,Dai,Fae,Fmi,Fme)
        td = makeT2(True,nelec,nspin,fs,spinints,ts,td,Dabij,Fae,Fmi,Fme,Wmnij,Wabef,Wmbej)
        ECCSD = ccsdenergy(nelec,nspin,fs,spinints,ts,td)
        DECC = abs(ECCSD - OLDCC)
        log('{0:<20d} {1:<20.6f} {2:<20.6f}'.format(Iter, ECCSD, DECC))
    log.hline()
    log('CCSD Energy = {}'.format(ECCSD))
    log.hline()
    log('Perturbative triples correction to CCSD'.format())
    Et = E_triples_correction(nelec,nspin,fs,ts,td,spinints)
    log('CCSD Triples Correction Energy = {}'.format(Et))
    log.hline()

    results = {
    "success":True,
    "CCSD_energy":ECCSD,
    "CCSD_T_energy":Et,
    }

    return results
