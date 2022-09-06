from moha.io.log import log,timer

import numpy as np

def rhf(S,Hcore,Eri,ndocc,I_thld,E_thld,D_thld):
    """Restricted close-shell Hartree Fock algorithm.

    Parameters
    ----------
    S: np.2darray
        Overlap matrix.

    Hcore: np.2darray
        Core Hamiltonian matrix.

    Eri: np.4darray
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
    
    #Store initial value of electronic energy and density matrix
    Dold = D
    Eold = np.einsum('pq,pq->', 2*Hcore, D)

    #Iteration section
    log.hline()
    log('{0:<20s} {1:<20s} {2:<20s} {3:<20s}'.format('Iter', 'E(elec)', 'dE','dRMS'))
    log.hline()
    for Iter in range(1,I_thld+1):
        #Update Fock matrix
        J = np.einsum('pqrs,rs->pq', Eri, D)
        K = np.einsum('prqs,rs->pq', Eri, D)
        F = Hcore + 2*J - K

        #Update density matrix
        Fp = X.dot(F).dot(X)
        Emo, Cp = np.linalg.eigh(Fp)
        C = X.dot(Cp)
        D = np.einsum('pi,qi->pq', C[:,:ndocc], C[:,:ndocc])

        #Update electronic energy
        E = np.einsum('pq,pq->', F + Hcore, D)

        #Update convergence
        dE = abs(E-Eold)
        dRMS = np.mean((D-Dold)**2)**0.5

        #Store Updated electronic energy and density matrix
        Dold = D
        Eold = E

        log('{0:<20d} {1:<20.6f} {2:<20.6f} {3:<20.6f}'.format(Iter, E, dE, dRMS))
        if (dE < E_thld) and (dRMS < D_thld):
            break
        if Iter == I_thld:
            raise Exception("Maximum number of SCF cycles exceeded.")

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
