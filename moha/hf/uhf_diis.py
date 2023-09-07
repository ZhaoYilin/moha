from moha.io.log import log,timer

import numpy as np

def uhf_diis(S,Hcore,Eri,na,nb,I_thld,E_thld,D_thld):
    """Restricted Hartree Fock algorithm.

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
    
    # Store initial value of electronic energy and density matrix
    Dold = D
    Eold = np.einsum('pq,pq->', 2*Hcore, D)

    # Iteration section
    log.hline()
    log('{0:<20s} {1:<20s} {2:<20s} {3:<20s}'.format('Iter', 'E(elec)', 'dE','dRMS'))
    log.hline()
    Fa_list = []
    Fb_list = []
    DIISA_error = []
    DIISB_error = []
    for Iter in range(1,I_thld+1):
        #Update Fock matrix
        Ja = np.einsum('pqrs,rs->pq', Eri, Da)
        Jb = np.einsum('pqrs,rs->pq', Eri, Db)
        Ka = np.einsum('prqs,rs->pq', Eri, Da)
        Kb = np.einsum('prqs,rs->pq', Eri, Db)
        Fa = Hcore + Ja+Jb - Ka
        Fb = Hcore + Ja+Jb - Kb

        Fa_list.append(Fa)
        Fb_list.append(Fb)

        #Update electronic energy
        E = np.einsum('pq,pq->', Hcore, D)
        E += np.einsum('pq,pq->', Fa, Da)
        E += np.einsum('pq,pq->', Fb, Db)
        E *= 0.5

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
        diisa_e = np.einsum('ij,jk,kl->il',Fa,Da,S)
        diisa_e -= np.einsum('ij,jk,kl->il',S,Da,Fa)
        diisa_e = X.dot(diisa_e).dot(X)

        diisb_e = np.einsum('ij,jk,kl->il',Fb,Db,S)
        diisb_e -= np.einsum('ij,jk,kl->il',S,Db,Fb)
        diisb_e = X.dot(diisb_e).dot(X)

        DIISA_error.append(diisa_e)
        DIISB_error.append(diisb_e)
        if Iter >= 2:
            # Limit size of DIIS vector
            diis_count = len(Fa_list)
            if diis_count > 6:
                # Remove oldest vector
                del Fa_list[0]
                del Fb_list[0]
                del DIISA_error[0]
                del DIISB_error[0]
                diis_count -= 1

            # Build error matrix B
            Ba = np.empty((diis_count + 1, diis_count + 1))
            Ba[-1, :] = -1
            Ba[:, -1] = -1
            Ba[-1, -1] = 0
            for num1, e1 in enumerate(DIISA_error):
                for num2, e2 in enumerate(DIISA_error):
                    if num2 > num1: continue
                    val = np.einsum('ij,ij->', e1, e2)
                    Ba[num1, num2] = val
                    Ba[num2, num1] = val

            Bb = np.empty((diis_count + 1, diis_count + 1))
            Bb[-1, :] = -1
            Bb[:, -1] = -1
            Bb[-1, -1] = 0
            for num1, e1 in enumerate(DIISB_error):
                for num2, e2 in enumerate(DIISB_error):
                    if num2 > num1: continue
                    val = np.einsum('ij,ij->', e1, e2)
                    Bb[num1, num2] = val
                    Bb[num2, num1] = val

            # normalize
            Ba[:-1, :-1] /= np.abs(Ba[:-1, :-1]).max()
            Bb[:-1, :-1] /= np.abs(Bb[:-1, :-1]).max()

            # Build residual vector, [Pulay:1980:393], Eqn. 6, RHS
            resid = np.zeros(diis_count + 1)
            resid[-1] = -1

            # Solve Pulay equations, [Pulay:1980:393], Eqn. 6
            cia = np.linalg.solve(Ba, resid)
            cib = np.linalg.solve(Bb, resid)

            # Calculate new fock matrix as linear
            # combination of previous fock matrices
            Fa = np.zeros_like(Fa)
            Fb = np.zeros_like(Fb)
            for num, c in enumerate(cia[:-1]):
                Fa += c * Fa_list[num]
            for num, c in enumerate(cib[:-1]):
                Fb += c * Fb_list[num]

        #Update density matrix
        Fp = X.dot(Fa).dot(X)
        Emo_a, Cp = np.linalg.eigh(Fp)
        Ca = X.dot(Cp)
        Da = np.einsum('pi,qi->pq', Ca[:,0:na], Ca[:,0:na])

        Fp = X.dot(Fb).dot(X)
        Emo_b, Cp = np.linalg.eigh(Fp)
        Cb = X.dot(Cp)
        Db = np.einsum('pi,qi->pq', Cb[:,0:nb], Cb[:,0:nb])


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
