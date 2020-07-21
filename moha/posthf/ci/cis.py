from moha.posthf.ci.auxiliary import CI_Basis_Set, CI_Hamiltonian
from moha.log import log,timer
import numpy as np
import copy

__all__ = ['CISSolver']

class CISSolver(object):
    """Configuration Interaction Singles Wavefunction.
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
        log('CIS Calculation'.format())
        log.hline()
        ham = copy.deepcopy(self.ham)
        hfwavefunction = copy.deepcopy(self.hfwavefunction)
        ci_basis = CI_Basis_Set(hfwavefunction,[0,1])
        ci_basis.generate_basis_set([0,1])
        H = CI_Hamiltonian(ham,hfwavefunction,ci_basis)
        ci_matrix = H.generate_matrix(ci_basis)
        e_cis, vec_cis = np.linalg.eigh(ci_matrix)
        E_cis_elec = e_cis[0]
        E_hf_elec = hfwavefunction.Eelec
        E_corr = E_cis_elec - E_hf_elec
        E_cis_tot = E_cis_elec + ham.operators['nuclear_repulsion'].value
        log('E_CIS = {}'.format(E_cis_tot))
        log('CCS Correlation energy = {}'.format(E_corr))
        log.hline()
        log('CIS Excitation Energies (Singlets only):')
        log.hline()
        log('#        Hartree')
        for i in range(1, len(e_cis)):
            excit_e = e_cis[i] + ham.operators['nuclear_repulsion'].value - hfwavefunction.Eelec
            log('{}      {}'.format(i, excit_e))
        return E_cis_tot

