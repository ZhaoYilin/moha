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
    Fprime = np.dot(np.transpose(X),np.dot(F,X))
    return Fprime


def C_matrix(F,X,Norb):
    """
    orbital energies and Coefficients matrix
    """
    Fprime = F_transform(F,X,Norb)
    eorbital,Cprime = np.linalg.eigh(Fprime)  
    C = np.dot(X,Cprime)  
    return eorbital,C


def D_matrix(C,Norb,occ):
    """
    Density matrix
    """
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
    F = np.zeros((Norb,Norb))
    for i in range(Norb):
        for j in range(Norb):
            F[i][j] = Hcore[i][j]+G[i][j]
    return F

# Calculate change in density matrix
def deltap(D,OLDD,nspatial):
    DELTA = 0.0e0
    for i in range(0,nspatial):
        for j in range(0,nspatial):
            DELTA = DELTA+((D[i,j]-OLDD[i,j])**2)
    DELTA = (DELTA/4)**(0.5)
    return DELTA

def energy(D,Hcore,F,Norb):  
    E = 0.
    for i in range(Norb):
        for j in range(Norb):
            E = E + 0.5*D[i][j]*(Hcore[i][j] + F[i][j])
    return E
