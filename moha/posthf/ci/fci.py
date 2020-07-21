from moha.posthf.ci.auxiliary import CI_Basis_Set, CI_Hamiltonian
from moha.log import log,timer
import numpy as np
import copy

class FullCISolver(object):
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
        log('Full CI Calculation'.format())
        log.hline()
        ham = copy.deepcopy(self.ham)
        hfwavefunction = copy.deepcopy(self.hfwavefunction)
        dim = hfwavefunction.dim *2
        occ = hfwavefunction.occ['alpha'] + hfwavefunction.occ['beta']
        vir = dim - occ
        excitation_level = list(range(min(occ,vir)+1))
        ci_basis = CI_Basis_Set(hfwavefunction,excitation_level)
        ci_basis.generate_basis_set(excitation_level)
        H = CI_Hamiltonian(ham,hfwavefunction,ci_basis)
        ci_matrix = H.generate_matrix(ci_basis)
        e_fci, vec_fci = np.linalg.eigh(ci_matrix)
        E_fci_elec = e_fci[0]
        E_hf_elec = hfwavefunction.Eelec
        E_corr = E_fci_elec - E_hf_elec
        E_fci_tot = E_fci_elec + ham.operators['nuclear_repulsion'].value
        log('E_Full_CI = {}'.format(E_fci_tot))
        log('E_corr = {}'.format(E_corr))
        return E_fci_tot

