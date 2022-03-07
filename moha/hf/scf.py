from moha.hf.auxiliary import *
from moha.hf.hf_wavefunction import HFWaveFunction
from moha.system.operator.base import OperatorNames 
from moha.io.log import log,timer

import copy

__all__ = ['PlainSCFSolver']

class PlainSCFSolver(object):
    """Plain self-consistent mean field method solver.

    Attributes
    ----------
    ham
        Chemical Hamiltonian.

    wfn
        Hartree Fock wavefunction.

    maxiter : int
        Maximum numer of iteration.

    E_conv : float
        Energy for convergence.

    D_conv : float
        Density for convergence.

    Methods
    -------
    __init__(self,ham,wfn,maxiter=50,E_conv=1.0E-8,D_conv=1.0E-3)
        Initialize the solver.

    kernel(self)
        Kernel of the solver.
    
    assign_hamiltonian(self,ham)
        Assign the chemical Hamiltonian to the solver.
    
    assign_wavefunction(self,wfn)
        Assign the Hartree Fock wavefunction to the solver.
    
    assign_maximum_iteration(self,maxiter)
        Assign the maximum number of iteration to the solver.
    
    assign_energy_convergence(self,E_conv)
        Assign the energy for convergence to the solver.
    
    assign_dentsity_convergence(self,D_conv)
        Assign the density for convergence to the solver.
    """
    def __init__(self,ham,wfn,maxiter=50,E_conv=1.0E-8,D_conv=1.0E-3):
        """Initialize the solver.

        Attributes
        ----------
        ham
            Chemical Hamiltonian.

        wfn
            Hartree Fock wavefunction.

        maxiter : int
            Maximum numer of iteration.

        E_conv : float
            Energy for convergence.

        D_conv : float
            Density for convergence.
        """
        self.assign_hamiltonian(ham)
        self.assign_wavefunction(wfn)
        self.assign_maximum_iteration(maxiter)
        self.assign_energy_convergence(E_conv)
        self.assign_density_convergence(D_conv)    
    
    @timer.with_section('SCF')
    def kernel(self):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            Hartree Fock calculation results.
        """
        log.hline()
        log('SCF Calculation Section'.format())
        log.hline()

        # Set defaults
        ham = self.ham
        occ = self.wfn.occ
        Norb = ham.nspatial
        
        maxiter = self.maxiter
        E_conv = self.E_conv
        D_conv = self.D_conv

        #Form the core Hamiltonian
        Hcore = ham.operators[OperatorNames.Hcore]
        #Build the orthogonalization matrix
        X = orthogonalization_matrix(ham.operators[OperatorNames.S])
        #Initial guess of density matrix
        D = np.zeros((Norb,Norb))

        #iteration section
        Enuc = ham.operators[OperatorNames.Enuc]
        RMS = 1.0
        Eold = 0.0
        log.hline()
        log('{0:2s}\t {1:3s}\t {2:4s}\t\t {3:5s}'.format('Iter', 'E(elec)', 'dE','dRMS'))
        log.hline()
        for Iter in range(1,maxiter+1):
            Iter += 1
            G = G_matrix(D,ham.operators[OperatorNames.Eri],Norb)
            F = F_matrix(Hcore,G,Norb)
            Eorbs,C = C_matrix(F,X,Norb)
            OLDD = D
            
            D = D_matrix(C,Norb,occ)
            dRMS = deltap(D,OLDD,Norb)

            Eelec = energy(D,Hcore,F,Norb)
            log('{0:2d}\t {1:3f}\t {2:4f}\t {3:5f}'.format(Iter, Eelec, Eelec-Eold, dRMS))
            if (abs(Eelec - Eold) < E_conv) and (dRMS < D_conv):
                break
            Eold = Eelec
            if Iter == maxiter:
                raise Exception("Maximum number of SCF cycles exceeded.")
        
        self.wfn.assign_coefficients(C)
        self.wfn.assign_density_matrix(D)
        self.wfn.assign_orbital_energies(Eorbs)
        
        log.hline()
        log('SCF Results'.format())
        log.hline()
        log('Electronic Energy = {}'.format(Eelec))
        log('Nuclear Energy = {}'.format(Enuc))
        log('Total Energy = {}'.format(Eelec+Enuc))
        log.hline()
        
        results = {
        "success": True,
        "electronic_energy": Eelec,
        "nuclear_energy": Enuc,
        "total_energy": Eelec+Enuc
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

    def assign_density_convergence(self,D_conv):
        """Assign the density for convergence to the solver.

        Attributes
        ----------
        D_conv : float
            Density for convergence.
        
        Raises
        ------
        TypeError
            If Density for convergence is not an integer.

        ValueError
            If Density for convergence is not a positive number.
        """
        if not isinstance(D_conv, float):
            raise TypeError("Density for convergence must be a float")
        if D_conv <= 0:
            raise ValueError("Density for convergence must be a positive number")        
        self.D_conv = D_conv
