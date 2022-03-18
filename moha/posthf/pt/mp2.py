from moha.system.hamiltonian.chemical_hamiltonian import *
from moha.system.operator.base import OperatorNames
from moha.posthf.pt.auxiliary import *
from moha.io.log import log, timer

import numpy as np
import copy

__all__ = ['MP2Solver']

class MP2Solver(object):
    """Second-order Moller-Plesset perturbation solver.
     
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

    @timer.with_section('MP2')
    def kernel(self):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            MP2 calculation results.    
        """
        log.hline()
        log('MP2 calculation section'.format())
        log.hline()

        ham = copy.deepcopy(self.ham)
        wfn = copy.deepcopy(self.wfn)
        hf_results = self.hf_results
        
        nspatial = ham.nspatial
        occ = wfn.occ
        C = wfn.coefficients
        eorbitals = wfn.orbital_energies

        Emp2 = 0.0
        if isinstance(ham,ChemicalHamiltonian):
            Eri = ham.operators[OperatorNames.Eri]
            Eri.basis_transformation(C)
            for i in range(occ['alpha']):
                for j in range(occ['alpha']):
                    for a in range(occ['alpha'],nspatial):
                        for b in range(occ['alpha'],nspatial):
                            Emp2 += Eri[i,a,j,b]*(2*Eri[i,a,j,b]-Eri[i,b,j,a])\
                            /(eorbitals[i] + eorbitals[j] -eorbitals[a] - eorbitals[b])
        
        log.hline()
        log('MP2 Results'.format())
        log.hline()
        log('{0:2s} {1:3f}'.format('Escf', hf_results['total_energy']))
        log('{0:2s} {1:3f}'.format('Emp2', Emp2))
        log('{0:2s} {1:3f}'.format('Etot', hf_results['total_energy']+Emp2))
        log.hline()
        
        results = {
        "success": True,
        "mp2_energy":Emp2,
        "total_energy":hf_results['total_energy']+Emp2
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
