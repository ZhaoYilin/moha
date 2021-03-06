from moha.posthf.ci.ci_basis_set import CIBasisSet
from moha.posthf.ci.ci_hamiltonian import RestrictedCIHamiltonian
from moha.log import log,timer
import numpy as np
import copy

__all__ = ['CISSolver']

class CISSolver(object):
    """Configuration Interaction Singles Wavefunction.
       CI with HF Ground state and all single and double excitations

    Attributes
    ----------
    ham
        Chemical Hamiltonian.

    wfn
        Hartree Fock wavefunction.

    hf_results : dict
        Hartree Fock calculation results.

    Methods
    -------
    __init__(self,ham,wfn,hf_results)
        Initialize the solver.

    kernel(self)
        Kernel of the solver.

    assign_hamiltonian(self,ham)
        Assign the chemical Hamiltonian to the solver.

    assign_wavefunction(self,wfn)
        Assign the Hartree Fock wavefunction to the solver.

    assign_hartree_fock_results(self,hf_results)
        Assign the Hartree Fock calculation results to the solver.
    """
    def __init__(self,ham,wfn,hf_results):
        """Initialize the solver.

        Attributes
        ----------
        ham
            Chemical Hamiltonian.

        wfn
            Hartree Fock wavefunction.

        hf_results : dict
            Hartree Fock calculation results.
        """
        self.assign_hamiltonian(ham)
        self.assign_wavefunction(wfn)
        self.assign_hartree_fock_results(hf_results)        
    
    @timer.with_section('CIS')
    def kernel(self):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            CIS calculation results.
        """
        log.hline()
        log('CIS Section'.format())
        log.hline()
        ham = copy.deepcopy(self.ham)
        wfn = copy.deepcopy(self.wfn)
        hf_results = self.hf_results

        ci_basis = CIBasisSet(wfn,[0,1])
        H = RestrictedCIHamiltonian(ham,wfn,ci_basis)
        ci_matrix = H.generate_matrix(ci_basis)
        e_ci, vec_ci = np.linalg.eigh(ci_matrix)
        E_ci_elec = e_ci[0]
        E_hf_elec = hf_results['electronic_energy']
        E_corr = E_ci_elec - E_hf_elec
        E_ci_tot = E_ci_elec + ham.operators['nuclear_repulsion'].integral
        
        log.hline()
        log('CIS Section'.format())
        log.hline()
        log('CIS Correlation energy = {}'.format(E_corr))
        log('Total energy = {}'.format(E_ci_tot))
        log.hline()
        
        results = {
        "total_energy":E_ci_tot
        }
        return results

    def assign_hamiltonian(self,ham):
        """Assign the chemical Hamiltonian to the solver.

        Attributes
        ----------
        ham
            Chemical Hamiltonian.
        """
        self.ham = ham

    def assign_wavefunction(self,wfn):
        """Assign the Hartree Fock wavefunction to the solver.

        Attributes
        ----------
        wfn
            Hartree Fock wavefunction.
        """
        self.wfn = wfn
    def assign_hartree_fock_results(self,hf_results):
        """Assign the Hartree Fock calculation results to the solver.

        Attributes
        ----------
        hf_results : dict
            Hartree Fock calculation results.

        Raises
        ------
        TypeError
            If Hartree Fock calculation results is not a dictionary.
        """
        if not isinstance(hf_results, dict):
            raise TypeError("Hartree Fock calculation results must be a dictionary")
        self.hf_results = hf_results
