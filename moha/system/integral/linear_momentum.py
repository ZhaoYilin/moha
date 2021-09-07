import numpy as np
from moha.system.integral.differential import D
from moha.system.integral.nuclear_attraction import E
def Pxyz(a,lmn1,A,b,lmn2,B):
    ''' Evaluates kinetic energy integral between two Gaussians
        Returns a float.
        a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
        lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
              for Gaussian 'a'
        lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
        A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
        B:    list containing origin of Gaussian 'b'
    '''
    i,k,m = lmn1
    j,l,n = lmn2
    term0 = D(i,j,1,A[0],B[0],a,b)*E(k,l,0,A[1],B[1],a,b)*E(m,n,0,A[2],B[2],a,b)
    term1 = E(i,j,0,A[0],B[0],a,b)*D(k,l,1,A[1],B[1],a,b)*E(m,n,0,A[2],B[2],a,b)
    term2 = E(i,j,0,A[0],B[0],a,b)*E(k,l,0,A[1],B[1],a,b)*D(m,n,1,A[2],B[2],a,b)
    return -1j*np.power(np.pi/(a+b),1.5)*(term0+term1+term2)

def linear_momentum(a,b):
    '''Evaluates kinetic energy between two contracted Gaussians
       Returns float.
       Arguments:
       a: contracted Gaussian 'a', BasisFunction object
       b: contracted Gaussian 'b', BasisFunction object
    '''
    t = 0.0
    for ia, ca in enumerate(a.coefs):
        for ib, cb in enumerate(b.coefs):
            t += a.norm[ia]*b.norm[ib]*ca*cb*\
                     Pxyz(a.exps[ia],a.shell,a.origin,\
                     b.exps[ib],b.shell,b.origin)
    return t

