from moha.system.operator.base import OperatorNames
from moha.posthf.ci.slater import SlaterDeterminant
from moha.posthf.ci.ci_basis_set import CIBasisSet
from moha.posthf.ci.ci_operator import CIOperator
from moha.posthf.ci.ci_wavefunction import CIWaveFunction
from moha.io.log import log,timer

import numpy as np
import copy

__all__ = ['FullCISolver']

class FullCISolver(object):
    """Full Configuration Interaction Wavefunction.
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
    
    @timer.with_section('FullCI')
    def kernel(self):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            Full CI calculation results.
        """
        log.hline()
        log('Full CI Section'.format())
        log.hline()
        ham = copy.deepcopy(self.ham)
        wfn = copy.deepcopy(self.wfn)
        hf_results = self.hf_results

        nelec = wfn.nelec
        nspatial = wfn.nspatial
        nspin = wfn.nspin
        

        nvir = nspin - nelec        
        excitation_level = min(nelec,nvir)+1

        reference = SlaterDeterminant.ground(nelec,nspin)
        cibs = CIBasisSet.build(reference,excitation_level)

        C = wfn.coefficients
        Hcore = ham.operators[OperatorNames.Hcore].basis_transformation(C)
        h1e = Hcore.spin_orbital_basis_integral
        Eri = ham.operators[OperatorNames.Eri].basis_transformation(C)
        g2e = Eri.double_bar

        ci_matrix = CIOperator.build(h1e, g2e, cibs)

        e_ci, vec_ci = np.linalg.eigh(ci_matrix)
        E_ci_elec = e_ci[0]
        E_hf_elec = hf_results['electronic_energy']
        E_corr = E_ci_elec - E_hf_elec
        E_ci_tot = E_ci_elec + ham.operators[OperatorNames.Enuc]

        #Build the ci wavefunction.
        ci_wfn = CIWaveFunction(nelec,nspatial,{},cibs,vec_ci[0])
        
        log.hline()
        log('Results'.format())
        log.hline()
        log('Full CI Correlation energy = {}'.format(E_corr))
        log('Total energy = {}'.format(E_ci_tot))
        log.hline()
        
        results = {
        "success":True,
        "CI_energy":E_ci_elec,
        "total_energy":E_ci_tot
        }

        return results, ci_wfn

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
