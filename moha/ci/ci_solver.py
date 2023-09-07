from moha.system.wavefunction.ci_wavefunction import CIWaveFunction
from moha.system.basis.basis import SlaterDeterminant
from moha.system.basis.basis_set import FockSpace
from moha.io.log import log,timer

import numpy as np

__all__ = ['CISolver']

class CISolver(object):
    """Configuration Interaction Singles and Doubles Wavefunction.
       CI with HF Ground state and all single and double excitations

    Attributes
    ----------
    ham
        Chemical Hamiltonian.

    wfn
        Hartree Fock wavefunction.

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
    """
    def __init__(self,ham,wfn):
        """Initialize the solver.

        Parameters
        ----------
        ham
            Chemical Hamiltonian.

        wfn
            Hartree Fock wavefunction.
        """
        self.assign_hamiltonian(ham)
        self.assign_wavefunction(wfn)

    @timer.with_section('CISD')
    def cis(self):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            CISD calculation results.
        """
        from moha.ci.cis import cis
        ham = self.ham
        wfn = self.wfn

        results = cis(ham,wfn)
        return results

    @timer.with_section('CISD')
    def cisd(self):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            CISD calculation results.
        """
        from moha.ci.cisd import cisd
        ham = self.ham
        wfn = self.wfn

        results = cisd(ham,wfn)
        return results
    
    @timer.with_section('CISD')
    def fci(self):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            CISD calculation results.
        """
        from moha.ci.fci import fci
        ham = self.ham
        wfn = self.wfn

        results = fci(ham,wfn)
        return results

    def assign_hamiltonian(self,ham):
        """Assign the chemical Hamiltonian to the solver.
        
        Parameters
        ----------
        ham: Hamiltonian
            Chemical Hamiltonian.
        """
        from moha.system.operator.hamiltonian import Hamiltonian
        if not isinstance(ham, Hamiltonian):
            raise TypeError("Parameter ham must be instance of Hamiltonian.")
        self.ham = ham

    def assign_wavefunction(self,wfn):
        """Assign the Hartree Fock wavefunction to the solver.
        
        Parameters
        ----------
        wfn: HFWaveFunction
            Hartree Fock wavefunction.
        """
        from moha.system.wavefunction.hf_wavefunction import HFWaveFunction
        if not isinstance(wfn, HFWaveFunction):
            raise TypeError("Parameter wfn must be instance of HFWaveFuntion.")
        self.wfn = wfn
