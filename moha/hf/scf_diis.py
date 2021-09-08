from moha.hf.auxiliary import *
from moha.system.wavefunction.hf_wavefunction import HFWaveFunction
from moha.io.log import log,timer

import copy

__all__ = ['DIISSCFSolver']

class DIISSCFSolver(object):
    """Direct Inversion in the Iterative Subspace (DIIS)" approach 
    self-consistent mean field method solver.

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
        log('SCF DIIS Calculation Section'.format())
        log.hline()
        # Set defaults
        ham = self.ham
        wfn = self.wfn
        maxiter = self.maxiter
        E_conv = self.E_conv
        D_conv = self.D_conv
        # Overlap Integral
        S = ham.operators['overlap'].integral

        # Get occ nbf and ndocc for closed shell molecules
        occ = wfn.occ
        nspatial = ham.nspatial
        ndocc = wfn.occ['alpha']

        # Compute required quantities for SCF
        V = ham.operators['nuclear_attraction'].integral
        T = ham.operators['kinetic'].integral
        I = ham.operators['electron_repulsion'].integral
        # Build H_core
        H = T + V

        # Orthogonalizer A = S^(-1/2)
        A = np.array(ham.operators['overlap'].integral)
        A = orthogonalization_matrix(A,type='symmetric')

        # Calculate initial core guess: [Szabo:1996] pp. 145
        Hp = A.dot(H).dot(A)            # Eqn. 3.177
        e, C2 = np.linalg.eigh(Hp)      # Solving Eqn. 1.178
        C = A.dot(C2)                   # Back transform, Eqn. 3.174
        Cocc = C[:, :ndocc]

        D = D_matrix(C,nspatial,occ)
        E = 0.0
        Enuc = ham.operators['nuclear_repulsion'].integral
        Eold = 0.0

        Fock_list = []
        DIIS_error = []
        log.hline()
        log('{0:2s}\t {1:3s}\t {2:4s}\t\t {3:5s}'.format('Iter', 'E(elec)', 'dE','dRMS'))
        log.hline()
        for Iter in range(1, maxiter + 1):

            # Build fock matrix
            J = np.einsum('pqrs,rs->pq', I, D)
            K = np.einsum('prqs,rs->pq', I, D)
            F = H + J * 2 - K

            # DIIS error build w/ HF analytic gradient ([Pulay:1969:197])
            diis_e = np.einsum('ij,jk,kl->il', F, D, S) - np.einsum('ij,jk,kl->il', S, D, F)
            diis_e = A.dot(diis_e).dot(A)
            Fock_list.append(F)
            DIIS_error.append(diis_e)
            dRMS = np.mean(diis_e**2)**0.5

            # SCF energy and update: [Szabo:1996], Eqn. 3.184, pp. 150
            Eelec = np.einsum('pq,pq->', F + H, D)
            log('{0:2d}\t {1:3f}\t {2:4f}\t {3:5f}'.format(Iter, Eelec, Eelec-Eold, dRMS))
            if (abs(Eelec - Eold) < E_conv) and (dRMS < D_conv):
                break

            Eold = Eelec

            if Iter >= 2:

                # Limit size of DIIS vector
                diis_count = len(Fock_list)
                if diis_count > 6:
                    # Remove oldest vector
                    del Fock_list[0]
                    del DIIS_error[0]
                    diis_count -= 1

                # Build error matrix B, [Pulay:1980:393], Eqn. 6, LHS
                B = np.empty((diis_count + 1, diis_count + 1))
                B[-1, :] = -1
                B[:, -1] = -1
                B[-1, -1] = 0
                for num1, e1 in enumerate(DIIS_error):
                    for num2, e2 in enumerate(DIIS_error):
                        if num2 > num1: continue
                        val = np.einsum('ij,ij->', e1, e2)
                        B[num1, num2] = val
                        B[num2, num1] = val

                # normalize
                B[:-1, :-1] /= np.abs(B[:-1, :-1]).max()

                # Build residual vector, [Pulay:1980:393], Eqn. 6, RHS
                resid = np.zeros(diis_count + 1)
                resid[-1] = -1

                # Solve Pulay equations, [Pulay:1980:393], Eqn. 6
                ci = np.linalg.solve(B, resid)

                # Calculate new fock matrix as linear
                # combination of previous fock matrices
                F = np.zeros_like(F)
                for num, c in enumerate(ci[:-1]):
                    F += c * Fock_list[num]

            # Diagonalize Fock matrix
            Fp = A.dot(F).dot(A)
            Eorbs, C2 = np.linalg.eigh(Fp)
            C = A.dot(C2)
            Cocc = C[:, :ndocc]
            D = np.einsum('pi,qi->pq', Cocc, Cocc)
            if Iter == maxiter:
                raise Exception("Maximum number of SCF cycles exceeded.")

        log.hline()
        log('SCF Results'.format())
        log.hline()
        log('Electronic Energy = {}'.format(Eelec))
        log('Nuclear Energy = {}'.format(Enuc))
        log('Total Energy = {}'.format(Eelec+Enuc))
        log.hline()

        self.wfn.assign_coefficients(C)

        results = {
        "success": True,
        "electronic_energy": Eelec,
        "nuclear_energy": Enuc,
        "total_energy": Eelec+Enuc,
        "orbital_energies": Eorbs
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
