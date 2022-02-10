from moha.system.hamiltonian.chemical_hamiltonian import *
from moha.posthf.pt.auxiliary import *
from moha.io.log import log, timer

import numpy as np
import copy

__all__ = ['MP3Solver']

class MP3Solver(object):
    """Third-order Moller-Plesset perturbation solver.

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

    @timer.with_section("MP3")
    def kernel(self):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            MP3 calculation results.
        """    
        log.hline()
        log('MP3 Calculation Section'.format())
        log.hline()
        
        ham = copy.deepcopy(self.ham)
        wfn = copy.deepcopy(self.wfn)
        hf_results = self.hf_results
        
        nspatial = ham.nspatial 
        occ = wfn.occ
        C = wfn.coefficients
        eorbitals = wfn.orbital_energies
        
        Emp2 = 0.0
        ham.operators['electron_repulsion'].basis_transformation(C)
        Eri = ham.operators['electron_repulsion'].integral
        for i in range(occ['alpha']):
            for j in range(occ['alpha']):
                for a in range(occ['alpha'],nspatial):
                    for b in range(occ['alpha'],nspatial):
                        Emp2 += Eri[i,a,j,b]*(2*Eri[i,a,j,b]-Eri[i,b,j,a])\
                        /(eorbitals[i] + eorbitals[j] -eorbitals[a] - eorbitals[b])

        Emp3 = 0.0
        Eri = ham.operators['electron_repulsion'].double_bar
        for i in range(occ['alpha']):
            for j in range(occ['alpha']):
                for k in range(occ['alpha']):
                    for l in range(occ['alpha']):
                        for a in range(occ['alpha'],nspatial):
                            for b in range(occ['alpha'],nspatial):
                                Emp3 += (1/8.0)*Eri[i,j,a,b]*Eri[k,l,i,j]*Eri[a,b,k,l]\
                                /((eorbitals[i] + eorbitals[j] -eorbitals[a] - eorbitals[b])\
                                *(eorbitals[k] + eorbitals[l] -eorbitals[a] - eorbitals[b]))
        for i in range(occ['alpha']):
            for j in range(occ['alpha']):
                for a in range(occ['alpha'],nspatial):
                    for b in range(occ['alpha'],nspatial):
                        for c in range(occ['alpha'],nspatial):
                            for d in range(occ['alpha'],nspatial):
                                Emp3 += (1/8.0)*Eri[i,j,a,b]*Eri[a,b,c,d]*Eri[c,d,i,j]\
                                /((eorbitals[i] + eorbitals[j] -eorbitals[a] - eorbitals[b])\
                                *(eorbitals[i] + eorbitals[j] -eorbitals[c] - eorbitals[d]))
        for i in range(occ['alpha']):
            for j in range(occ['alpha']):
                for k in range(occ['alpha']):
                    for a in range(occ['alpha'],nspatial):
                        for b in range(occ['alpha'],nspatial):
                            for c in range(occ['alpha'],nspatial):
                                Emp3 += Eri[i,j,a,b]*Eri[k,b,c,j]*Eri[a,c,i,k]\
                                /((eorbitals[i] + eorbitals[j] -eorbitals[a] - eorbitals[b])\
                                *(eorbitals[i] + eorbitals[k] -eorbitals[c] - eorbitals[c]))
        
        log.hline()
        log('MP3 Results'.format())
        log.hline()
        log('{0:2s} {1:3f}'.format('Escf', hf_results['total_energy']))
        log('{0:2s} {1:3f}'.format('Emp2', Emp2))
        log('{0:2s} {1:3f}'.format('Emp3', Emp3))
        log('{0:2s} {1:3f}'.format('Etot', hf_results['total_energy']+Emp2+Emp3))
        log.hline()

        results = {
        "success": True,
        "MP2_energy":Emp2,
        "MP3_energy":Emp3,
        "total_energy":hf_results['total_energy']+Emp2+Emp3
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
