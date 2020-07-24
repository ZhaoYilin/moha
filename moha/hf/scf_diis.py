from moha.hf.auxiliary import *
from moha.log import log,timer
from moha.system.state import HFWaveFunction

import copy

__all__ = ['SCFDIISSolver']

class SCFDIISSolver(object):
    """
    """
    def __init__(self,ham,occ,maxiter=100,cutoff=1e-6):
        """
        """
        self.ham = ham
        self.occ = occ
        self.maxiter = maxiter
        self.cutoff = cutoff

    def kernel(self):
        """
        """
        log.hline()
        log('SCF DIIS Calculation'.format())
        log.hline()
        # Set defaults
        ham = copy.deepcopy(self.ham)
        maxiter = 40
        E_conv = 1.0E-8
        D_conv = 1.0E-3

        # Integral generation from Psi4's MintsHelper
        S = ham.operators['overlap'].value

        # Get nbf and ndocc for closed shell molecules
        nbf = ham.dim
        ndocc = self.occ['alpha']

        print('\nNumber of occupied orbitals: %d' % ndocc)
        print('Number of basis functions: %d' % nbf)


        # Compute required quantities for SCF
        V = ham.operators['nuclear_attraction'].value
        T = ham.operators['kinetic'].value
        I = ham.operators['electron_repulsion'].tensor_format()
        # Build H_core
        H = T + V

        # Orthogonalizer A = S^(-1/2)
        A = np.array(ham.operators['overlap'].value)
        A = orthogonalization_matrix(A,type='symmetric')

        # Calculate initial core guess: [Szabo:1996] pp. 145
        Hp = A.dot(H).dot(A)            # Eqn. 3.177
        e, C2 = np.linalg.eigh(Hp)      # Solving Eqn. 1.178
        C = A.dot(C2)                   # Back transform, Eqn. 3.174
        Cocc = C[:, :ndocc]

        D = D_matrix(C,nbf,{'alpha':5,'beta':5})

        print('\nStart SCF iterations:\n')
        E = 0.0
        Enuc = ham.operators['nuclear_repulsion'].value
        Eold = 0.0

        Fock_list = []
        DIIS_error = []
        for SCF_ITER in range(1, maxiter + 1):

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
            SCF_E = np.einsum('pq,pq->', F + H, D) + Enuc

            print('SCF Iteration %3d: Energy = %4.16f   dE = % 1.5E   dRMS = %1.5E'
                % (SCF_ITER, SCF_E, (SCF_E - Eold), dRMS))
            if (abs(SCF_E - Eold) < E_conv) and (dRMS < D_conv):
                break

            Eold = SCF_E

            if SCF_ITER >= 2:

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
            e, C2 = np.linalg.eigh(Fp)
            C = A.dot(C2)
            Cocc = C[:, :ndocc]
            D = np.einsum('pi,qi->pq', Cocc, Cocc)

            if SCF_ITER == maxiter:
                raise Exception("Maximum number of SCF cycles exceeded.")


        print('Final SCF energy: %.8f hartree' % SCF_E)

        return 0
        #return HFWaveFunction(Norb,occ,C,D,F,eorbitals,Eelec,Etot)

