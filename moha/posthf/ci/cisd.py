from moha.posthf.ci.auxiliary import CI_Basis_Set, CI_Hamiltonian
from moha.log import log,timer
import numpy as np
import copy

__all__ = ['CISDSolver']

class CISDSolver(object):
    """Configuration Interaction Singles and Doubles Wavefunction.
    CI with HF Ground state and all single and double excitations
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
        log('CISD'.format())
        log.hline()
        ham = copy.deepcopy(self.ham)
        hfwavefunction = copy.deepcopy(self.hfwavefunction)
        ci_basis = CI_Basis_Set(hfwavefunction,[0,1,2])
        ci_basis.generate_basis_set([0,1,2])
        H = CI_Hamiltonian(ham,hfwavefunction,ci_basis)
        ci_matrix = H.generate_matrix(ci_basis)
        e_cisd, vec_cisd = np.linalg.eigh(ci_matrix)
        E_cisd_elec = e_cisd[0]
        E_hf_elec = hfwavefunction.Eelec
        E_corr = E_cisd_elec - E_hf_elec
        E_cisd_tot = E_cisd_elec + ham.operators['nuclear_repulsion'].value
        log('E_CISD = {}'.format(E_cisd_tot))
        log('CCSD Correlation energy = {}'.format(E_corr))
        return E_cisd_tot

