from moha.posthf.cc.auxiliary import *
from moha.log import log, timer

import numpy as np
import copy

__all__ = ['CCSDSolver']

class CCSDSolver(object):
    """The coupled cluster singles and double solver.

    Attributes
    ----------
    ham
        Chemical Hamiltonian.

    wfn
        Hartree Fock wavefunction.

    hf_results : dict
        Hartree Fock calculation results.

    maxiter : int
        Maximum numer of iteration.

    E_conv : float
        Energy for convergence.            
    
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

    assign_maximum_iteration(self,maxiter)
        Assign the maximum number of iteration to the solver.

    assign_energy_convergence(self,E_conv)
        Assign the energy for convergence to the solver.        
    """
    def __init__(self,ham,wfn,hf_results,maxiter=100,E_conv=1.0E-8):
        """Initialize the solver.

        Attributes
        ----------
        ham
            Chemical Hamiltonian.

        wfn
            Hartree Fock wavefunction.

        hf_results : dict
            Hartree Fock calculation results.
        
        maxiter : int
            Maximum numer of iteration.

        E_conv : float
            Energy for convergence.            
        """
        self.assign_hamiltonian(ham)
        self.assign_wavefunction(wfn)
        self.assign_hartree_fock_results(hf_results)        
        self.assign_maximum_iteration(maxiter)
        self.assign_energy_convergence(E_conv)

    @timer.with_section('CCSD')
    def kernel(self):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            CCSD calculation results.
        """
        log.hline()
        log('CCSD Section'.format())
        log.hline()
        ham = copy.deepcopy(self.ham)
        wfn = copy.deepcopy(self.wfn)
        hf_results = self.hf_results

        occ = wfn.occ
        Nelec = occ['alpha'] + occ['beta']
        C = wfn.coefficients
        nspin = ham.nspin
        #Transfer Fock integral from spatial to spin basis
        fs = spinfock(hf_results['orbital_energies'])
        #Transfer electron repulsion integral from atomic basis
        #to molecular basis
        ham.operators['electron_repulsion'].basis_transformation(C)
        #build double bar integral <ij||kl>
        spinints = ham.operators['electron_repulsion'].double_bar 
        ts,td,Dai,Dabij = initialize(fs,spinints,Nelec,nspin)
        ECCSD = 0
        DECC = 1.0
        Iter = 0
        log.hline()
        log('{0:2s} {1:3s} {2:4s}'.format('Iter', 'ECC', 'Delta'))
        log.hline()
        while DECC > self.E_conv or Iter> self.maxiter: # arbitrary convergence criteria
            Iter += 1
            OLDCC = ECCSD
            Fae,Fmi,Fme,Wmnij,Wabef,Wmbej = updateintermediates(True,Nelec,nspin,fs,spinints,ts,td)
            ts = makeT1(True,Nelec,nspin,fs,spinints,ts,td,Dai,Fae,Fmi,Fme)
            td = makeT2(True,Nelec,nspin,fs,spinints,ts,td,Dabij,Fae,Fmi,Fme,Wmnij,Wabef,Wmbej)
            ECCSD = ccsdenergy(Nelec,nspin,fs,spinints,ts,td)
            DECC = abs(ECCSD - OLDCC)
            log('{0:2d} {1:3f} {2:4f}'.format(Iter, ECCSD, DECC))
        log.hline()
        log('CCSD Energy = {}'.format(ECCSD))
        log('Total Energy = {}'.format(hf_results['total_energy']+ECCSD))        
        log.hline()

        results = {
        "success":True,
        "CCSD_energy":ECCSD,
        "total_energy":hf_results['total_energy']+ECCSD
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

    def assign_maximum_iteration(self,maxiter):
        """Assign the maximum number of iteration to the solver.

        Attributes
        ----------
        maxiter : int
            Maximum numer of iteration.

        Raises
        ------
        TypeError
            If maximum number of iteration is not an integer.

        ValueError
            If maximum number of iteration is not a positive number.
        """
        if not isinstance(maxiter, int):
            raise TypeError("Maximum number of iteration must be an integer")
        if maxiter <= 0:
            raise ValueError("Maximum number of iteration must be a positive number")
        self.maxiter = maxiter

    def assign_energy_convergence(self,E_conv):
        """Assign the energy for convergence to the solver.

        Attributes
        ----------
        E_conv : float
            Energy for convergence.

        Raises
        ------
        TypeError
            If Energy for convergence is not an integer.

        ValueError
            If Energy for convergence is not a positive number.
        """
        if not isinstance(E_conv, float):
            raise TypeError("Energy for convergence must be a float")
        if E_conv <= 0:
            raise ValueError("Energy for convergence must be a positive number")
        self.E_conv = E_conv        
