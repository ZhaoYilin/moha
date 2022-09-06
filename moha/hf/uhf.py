from moha.io.log import log,timer

import numpy as np

def uhf(S,Hcore,Eri,na,nb,I_thld,E_thld,D_thld):
    """Restricted open-shell Hartree Fock algorithm.

    Parameters
    ----------
    S: np.ndarray
        Overlap matrix.

    Hcore: np.ndarray
        Core Hamiltonian matrix.

    Eri: np.ndarray
        Electron repuslion matrix.

    na: int
        Number of spin alpha electrons.

    nb: int
        Number of spin beta electrons.

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
    Da = np.einsum('pi,qi->pq', C[:,0:na], C[:,0:na])
    Db = np.einsum('pi,qi->pq', C[:,0:nb], C[:,0:nb])
    D = Da + Db
 
    #Store initial value of electronic energy and density matrix
    Dold = D
    Eold = np.einsum('pq,pq->', 2*Hcore, D)

    #Iteration section
    log.hline()
    log('{0:<20s} {1:<20s} {2:<20s} {3:<20s}'.format('Iter', 'E(elec)', 'dE','dRMS'))
    log.hline()
    for Iter in range(1,I_thld+1):
        #Update Fock matrix
        Ja = np.einsum('pqrs,rs->pq', Eri, Da)
        Jb = np.einsum('pqrs,rs->pq', Eri, Db)
        Ka = np.einsum('prqs,rs->pq', Eri, Da)
        Kb = np.einsum('prqs,rs->pq', Eri, Db)
        Fa = Hcore + Ja+Jb - Ka
        Fb = Hcore + Ja+Jb - Kb

        #Update density matrix
        Fp = X.dot(Fa).dot(X)
        Emo_a, Cp = np.linalg.eigh(Fp)
        Ca = X.dot(Cp)
        Da = np.einsum('pi,qi->pq', Ca[:,0:na], Ca[:,0:na])

        Fp = X.dot(Fb).dot(X)
        Emo_b, Cp = np.linalg.eigh(Fp)
        Cb = X.dot(Cp)
        Db = np.einsum('pi,qi->pq', Cb[:,0:nb], Cb[:,0:nb])

        D = Da + Db

        #Update electronic energy
        E = np.einsum('pq,pq->', Hcore, D)
        E += np.einsum('pq,pq->', Fa, Da)
        E += np.einsum('pq,pq->', Fb, Db)
        E *= 0.5

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

    Emo = np.block([Emo_a,Emo_b])
    C = np.block([Ca,Cb])
    D = np.block([Da,Db])

    results = {
    "success": True,
    "energy": E,
    "Emo": Emo,
    "C": C,
    "D": D
    }
    return results
