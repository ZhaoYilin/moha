from moha.io.log import log,timer

import numpy as np

def rhf_diis(S,Hcore,Eri,ndocc,I_thld,E_thld,D_thld):
    """Restricted close-shell Hartree Fock algorithm.

    Parameters
    ----------
    S: np.ndarray
        Overlap matrix.

    Hcore: np.ndarray
        Core Hamiltonian matrix.

    Eri: np.ndarray
        Electron repuslion matrix.

    ndocc: int
        Number of doubly occupied orbitals.

    I_thld: int
        Threshold of iteration number.

    E_thld: float
        Threshold of energy convergence.

    D_thld: float
            Threshold of density matrix convergence.

    Returns
    -------
    results : dict
        Hartree Fock calculation results.
    """
    #Build the orthogonalization matrix
    val, vec = np.linalg.eig(S)  
    val_minus_half = (np.diag(val**(-0.5))) 
    X = np.dot(vec,np.dot(val_minus_half,np.transpose(vec)))

    #Initial guess of density matrix
    Fp = X.dot(Hcore).dot(X)
    Emo, Cp = np.linalg.eigh(Fp)
    C = X.dot(Cp)
    D = np.einsum('pi,qi->pq', C[:,:ndocc], C[:,:ndocc])
    
    # Store initial value of electronic energy and density matrix
    Dold = D
    Eold = np.einsum('pq,pq->', 2*Hcore, D)

    # Iteration section
    log.hline()
    log('{0:<20s} {1:<20s} {2:<20s} {3:<20s}'.format('Iter', 'E(elec)', 'dE','dRMS'))
    log.hline()
    F_list = []
    DIIS_error = []
    for Iter in range(1,I_thld+1):
        # Update Fock matrix
        J = np.einsum('pqrs,rs->pq', Eri, D)
        K = np.einsum('prqs,rs->pq', Eri, D)
        F = Hcore + 2*J - K
        F_list.append(F)

        # Update electronic energy
        E = np.einsum('pq,pq->', F + Hcore, D)

        # Update convergence
        dE = abs(E-Eold)
        dRMS = np.mean((D-Dold)**2)**0.5

        # Store Updated electronic energy and density matrix
        Dold = D
        Eold = E

        log('{0:<20d} {1:<20.6f} {2:<20.6f} {3:<20.6f}'.format(Iter, E, dE, dRMS))
        if dE < E_thld and dRMS < D_thld:
            break
        if Iter == I_thld:
            raise Exception("Maximum number of SCF cycles exceeded.")

        # DIIS section
        diis_e = np.einsum('ij,jk,kl->il',F,D,S)-np.einsum('ij,jk,kl->il',S,D,F)
        diis_e = X.dot(diis_e).dot(X)
        DIIS_error.append(diis_e)
        if Iter >= 2:
            # Limit size of DIIS vector
            diis_count = len(F_list)
            if diis_count > 6:
                # Remove oldest vector
                del F_list[0]
                del DIIS_error[0]
                diis_count -= 1

            # Build error matrix B
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
                F += c * F_list[num]

        # Update density matrix
        Fp = X.dot(F).dot(X)
        Emo, Cp = np.linalg.eigh(Fp)
        C = X.dot(Cp)
        D = np.einsum('pi,qi->pq', C[:,:ndocc], C[:,:ndocc])

    Emo = np.block([Emo,Emo])
    C = np.block([C,C])
    D = np.block([D,D])

    results = {
    "success": True,
    "energy": E,
    "Emo": Emo,
    "C": C,
    "D": D
    }
    return results
