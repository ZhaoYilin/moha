from moha.hf.auxiliary import *
from moha.log import log,timer
from moha.system.state import HFWaveFunction

import copy

__all__ = ['PlainSCFSolver']

class PlainSCFSolver(object):
    """
    """
    def __init__(self,ham,occ,maxiter=100,cutoff=1e-6):
        """
        """
        self.ham = ham
        self.occ = occ
        self.maxiter = maxiter
        self.cutoff = cutoff

    def kernel(self):
        """
        """
        log.hline()
        log('SCF Calculation'.format())
        log.hline()
        ham = copy.deepcopy(self.ham)
        occ = self.occ
        Norb = ham.dim
        #Form the core Hamiltonian
        Hcore = ham.operators['kinetic'].value + ham.operators['nuclear_attraction'].value
        #Build the orthogonalization matrix
        X = orthogonalization_matrix(ham.operators['overlap'].value)
        #Initial guess of density matrix
        if occ['alpha']==occ['beta']:
            D = np.zeros((Norb,Norb))
        else:
            D = {'alpha':np.zeros((Norb,Norb)),'beta':np.zeros((Norb,Norb))}
        #iteration section
        Iter = 0
        RMS = 1.0
        log.hline()
        log('{0:2s} {1:3s} {2:4s} {3:5s}'.format('Iter', 'E(elec)', 'E(tot)','RMS'))
        log.hline()
        while RMS > self.cutoff and Iter < self.maxiter:
            Iter += 1
            G = G_matrix(D,ham.operators['electron_repulsion'].value,Norb)
            F = F_matrix(Hcore,G,Norb)
            eorbitals,C = C_matrix(F,X,Norb)
            OLDD = D
            
            D = D_matrix(C,Norb,occ)
            RMS = deltap(D,OLDD,Norb)

            Eelec = energy(D,Hcore,F,Norb)
            Etot = Eelec + ham.operators['nuclear_repulsion'].value
            log('{0:2d} {1:3f} {2:4f} {3:5f}'.format(Iter, Eelec, Etot, RMS))
        log('Converged SCF energy = {}'.format(Etot))
        
        return HFWaveFunction(Norb,occ,C,D,F,eorbitals,Eelec,Etot)

