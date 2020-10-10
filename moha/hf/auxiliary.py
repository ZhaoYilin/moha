import math
import numpy as np  
####################################  
#                                  #
#            FUNCTIONS             #
#                                  #
####################################  

def eint(a,b,c,d):
    if a > b: ab = a*(a+1)/2 + b
    else: ab = b*(b+1)/2 + a
    if c > d: cd = c*(c+1)/2 + d
    else: cd = d*(d+1)/2 + c
    if ab > cd: abcd = ab*(ab+1)/2 + cd
    else: abcd = cd*(cd+1)/2 + ab
    return int(abcd)

"""
def tei(a,b,c,d):
    return twoe.get(eint(a,b,c,d),0.0)
"""

def orthogonalization_matrix(S,type='symmetric'):
    """
    Diagonlize the overlap matrix
    Schmidt
    symmetric
    canonical
    """
    if type == 'Schmidt':
        pass
    elif type == 'symmetric':
        val, vec = np.linalg.eig(S)  
        val_minus_half = (np.diag(val**(-0.5)))  
        X = np.dot(vec,np.dot(val_minus_half,np.transpose(vec)))  
    elif type == 'canonical':
        val, vec = np.linalg.eig(S)  
        val_minus_half = (np.diag(val**(-0.5)))  
        X = np.dot(vec,val_minus_half)  
    return X


def F_transform(F,X,Norb):
    if type(F) is dict:
        Fprime = {'alpha':np.zeros((Norb,Norb)),'beta':np.zeros((Norb,Norb))}
        for spin in F:
            Fprime[spin] = np.dot(np.transpose(X),np.dot(F[spin],X))
    elif type(F) is np.ndarray:
        Fprime = np.dot(np.transpose(X),np.dot(F,X))
    return Fprime


def C_matrix(F,X,Norb):
    """
    orbital energies and Coefficients matrix
    """
    if type(F) is dict:
        eorbital = {'alpha':np.zeros(Norb),'beta':np.zeros(Norb)}
        C = {'alpha':np.zeros((Norb,Norb)),'beta':np.zeros((Norb,Norb))}
        Cprime = {'alpha':np.zeros((Norb,Norb)),'beta':np.zeros((Norb,Norb))}
        Fprime = F_transform(F,X,Norb)
        for spin in F:
            eorbital[spin],Cprime[spin] = np.linalg.eigh(Fprime[spin]) 
            e,evec = np.linalg.eigh(Fprime[spin]) 
            eorbital[spin] = np.argsort(e.real)
            Cprime[spin] = evec[:,eorbital[spin]]

            C[spin] = np.dot(X,Cprime[spin])  
    elif type(F) is np.ndarray:
        Fprime = F_transform(F,X,Norb)
        eorbital,Cprime = np.linalg.eigh(Fprime)  
        C = np.dot(X,Cprime)  
    return eorbital,C


def D_matrix(C,Norb,occ):
    """
    Density matrix
    """
    if type(C) is dict:
        D = {'alpha':np.zeros((Norb,Norb)),'beta':np.zeros((Norb,Norb))}
        for spin in C:
            for i in range(Norb):
                for j in range(Norb):
                    for k in range(occ[spin]):
                        D[spin][i,j] += C[spin][i,k]*C[spin][j,k]
    elif type(C) is np.ndarray:
        D = np.zeros((Norb,Norb))
        for i in range(Norb):
            for j in range(Norb):
                for k in range(occ['alpha']):
                    D[i,j] += 2*C[i,k]*C[j,k]
    return D

def G_matrix(D,Eri,Norb):
    """
    G matrix
    """
    if type(D) is dict:
        Dtot = D['alpha']+D['beta']
        G = {'alpha':np.zeros((Norb,Norb)),'beta':np.zeros((Norb,Norb))}
        for spin in D:
            for i in range(Norb):
                for j in range(Norb):
                    for k in range(Norb):
                        for l in range(Norb):
                            G[spin][i,j] += Dtot[k,l]*Eri[i,j,k,l]-D[spin][k,l]*Eri[i,k,j,l]
    elif type(D) is np.ndarray:
        G = np.zeros((Norb,Norb))
        for i in range(Norb):
            for j in range(Norb):
                for k in range(Norb):
                    for l in range(Norb):
                        G[i][j] += D[k][l]*(Eri[i,j,k,l]-0.5*Eri[i,k,j,l])
    return G

def F_matrix(Hcore,G,Norb):
    """
    Fock matrix
    """
    if type(G) is dict:
        F = {'alpha':np.zeros((Norb,Norb)),'beta':np.zeros((Norb,Norb))}
        for spin in G:
            for i in range(Norb):
                for j in range(Norb):
                    F[spin][i,j] = Hcore[i,j]+G[spin][i,j]
    elif type(G) is np.ndarray:
        F = np.zeros((Norb,Norb))
        for i in range(Norb):
            for j in range(Norb):
                F[i][j] = Hcore[i][j]+G[i][j]
    return F

# Calculate change in density matrix
def deltap(D,OLDD,nspatial):
    DELTA = 0.0e0
    if type(D) is dict:
        D = D['alpha']+D['beta']
        OLDD = OLDD['alpha']+OLDD['beta']
    for i in range(0,nspatial):
        for j in range(0,nspatial):
            DELTA = DELTA+((D[i,j]-OLDD[i,j])**2)
    DELTA = (DELTA/4)**(0.5)
    return DELTA

def energy(D,Hcore,F,Norb):  
    E = 0.
    if type(D) is dict:
        Dtot = D['alpha']+D['beta']
        for i in range(Norb):  
            for j in range(Norb):  
                E = E + 0.5*(Dtot[i,j]*Hcore[i,j] + D['alpha'][i,j]*F['alpha'][i,j] + D['beta'][i,j]*F['beta'][i,j])  
    elif type(D) is np.ndarray:
        for i in range(Norb):
            for j in range(Norb):
                E = E + 0.5*D[i][j]*(Hcore[i][j] + F[i][j])
    return E
