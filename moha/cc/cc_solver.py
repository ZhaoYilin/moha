from moha.io.log import log, timer

import numpy as np
import copy

__all__ = ['CCSolver']

class CCSolver(object):
    """The coupled cluster singles and double solver.

    Attributes
    ----------
    ham: Hamiltonian
        Chemical Hamiltonian.

    wfn: HFWaveFunction
        Hartree Fock wavefunction.

    maxiter : int
        Maximum numer of iteration.

    E_conv : float
        Energy for convergence.            
    
    Methods
    -------
    __init__(self,ham,wfn)
        Initialize the solver.

    ccsd(self)
        Kernel of the solver.

    ccsd_t(self)
        Kernel of the solver.

    assign_hamiltonian(self,ham)
        Assign the chemical Hamiltonian to the solver.

    assign_wavefunction(self,wfn)
        Assign the Hartree Fock wavefunction to the solver.

    assign_maximum_iteration(self,maxiter)
        Assign the maximum number of iteration to the solver.

    assign_energy_convergence(self,E_conv)
        Assign the energy for convergence to the solver.        
    """
    def __init__(self,ham,wfn,I_thld=100,E_thld=1.0E-8):
        """Initialize the solver.

        Parameters
        ----------
        ham: Hamiltonian
            Chemical Hamiltonian.

        wfn: HFWaveFunction
            Hartree Fock wavefunction.

        maxiter : int
            Maximum numer of iteration.

        E_conv : float
            Energy for convergence.            
        """
        self.assign_hamiltonian(ham)
        self.assign_wavefunction(wfn)
        self.assign_iteration_threshold(I_thld)
        self.assign_energy_threshold(E_thld)

    @timer.with_section('CCSD')
    def ccsd(self):
        """Kernel of the solver.

        Return
        ------
        energy : float
            CC calculation results.
        """
        from moha.cc.ccsd import ccsd

        ham = self.ham
        wfn = self.wfn
        I_thld = self.I_thld
        E_thld = self.E_thld

        # Prepare parameters for mp2 algorithm
        Eri = ham['Eri']
        C = wfn.coefficient
        Emo = wfn.orbital_energies
        na = wfn.ansatz.nalphas
        nb = wfn.ansatz.nbetas

        # Split C and Emo to alpha and beta section
        Ca,Cb = np.hsplit(C,2)
        Emo_a,Emo_b = np.hsplit(Emo,2)


        E = ccsd(Eri,Ca,Cb,Emo_a,Emo_b,na,nb,I_thld,E_thld)

        return E

    def assign_hamiltonian(self,ham):
        """Assign the chemical Hamiltonian to the solver.
        
        Parameters
        ----------
        ham: Hamiltonian
            Chemical Hamiltonian.

        Raises
        ------
        TypeError
            If ham is not Hamiltonian instance.
        """
        from moha.hamiltonian.hamiltonian import Hamiltonian
        if not isinstance(ham, Hamiltonian):
            raise TypeError("Parameter ham must be Hamiltonian instance.")
        self.ham = ham

    def assign_wavefunction(self,wfn):
        """Assign the Hartree Fock wavefunction to the solver.
        
        Parameters
        ----------
        wfn: HFWaveFunction
            Hartree Fock wavefunction.
        """
        from moha.wavefunction.hf_wavefunction import HFWaveFunction
        if not isinstance(wfn, HFWaveFunction):
            raise TypeError("Parameter wfn must be instance of HFWaveFuntion.")
        self.wfn = wfn

    def assign_iteration_threshold(self,I_thld):
        """Assign the threshold of iteration number to the solver.

        Parameters
        ----------
        I_thld: int
            Threshold of iteration number.
        
        Raises
        ------
        TypeError
            If threshold of iteration number is not integer.

        ValueError
            If threshold of iteration number is not positive.
        """
        if not isinstance(I_thld, int):
            raise TypeError("Threshold of iteration number must be integer.")
        if I_thld <= 0:
            raise ValueError("Threshold of iteration number must be positive.")
        self.I_thld = I_thld

    def assign_energy_threshold(self,E_thld):
        """Assign the threshold of energy convergence to the solver.

        Parameters
        ----------
        E_thld: float
            Threshold of energy convergence.
        
        Raises
        ------
        TypeError
            If threshold of energy convergence is not integer.

        ValueError
            If threshold of energy convergence is not positive.
        """
        if not isinstance(E_thld, float):
            raise TypeError("Threshold of energy convergence must be float.")
        if E_thld <= 0:
            raise ValueError("Threshold of energy convergence must be positive.")
        self.E_thld = E_thld

