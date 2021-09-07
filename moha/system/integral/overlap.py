import numpy as np

def S(i,j,A,B,a,b):
    ''' Recursive definition of Hermite Gaussian coefficients.
        Returns a float.

    Parameters
    ----------
    a : 
        orbital exponent on Gaussian 'a' (e.g. alpha in the text)
    b : 
        orbital exponent on Gaussian 'b' (e.g. beta in the text)
    i,j: 
        orbital angular momentum number on Gaussian 'a' and 'b'
    t: 
        number nodes in Hermite (depends on type of integral,
           e.g. always zero for overlap integrals)
        Qx: distance between origins of Gaussian 'a' and 'b'
    '''
    p = a + b
    q = (a*b)/p
    P = (1./p)*(a*A + b*B)
    PA = P-A
    PB = P-B
    AB = A-B

    if i == j == 0:
        # base case
        return np.exp(-q*AB*AB) # K_AB

    elif i<0 or j<0 :
        return 0.

    elif i == 0 and j>0:
        # decrement index j
        return  PB*S(i,j-1,A,B,a,b) +\
                1./(2*p)*i*S(i-1,j-1,A,B,a,b) +\
                1./(2*p)*(j-1)*S(i,j-2,A,B,a,b)

    else:
        # decrement index i
        return  PA*S(i-1,j,A,B,a,b) +\
                1./(2*p)*(i-1)*S(i-2,j,A,B,a,b) +\
                1./(2*p)*j*S(i-1,j-1,A,B,a,b)

def Sxyz(a,lmn1,A,b,lmn2,B):
    ''' Evaluates overlap integral between two Gaussians
        Returns a float.

    Parameters
    ----------
        a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
        lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
              for Gaussian 'a'
        lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
        A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
        B:    list containing origin of Gaussian 'b'
    '''
    i,k,m = lmn1 # shell angular momentum on Gaussian 'a'
    j,l,n = lmn2 # shell angular momentum on Gaussian 'b'
    Sx = S(i,j,A[0],B[0],a,b) # X
    Sy = S(k,l,A[1],B[1],a,b) # Y
    Sz = S(m,n,A[2],B[2],a,b) # Z
    return Sx*Sy*Sz*np.power(np.pi/(a+b),1.5)

def overlap(a,b):
    '''Evaluates overlap between two contracted Gaussians
       Returns float.
       Arguments:
       a: contracted Gaussian 'a', BasisFunction object
       b: contracted Gaussian 'b', BasisFunction object
    '''
    s = 0.0
    for ia, ca in enumerate(a.coefs):
        for ib, cb in enumerate(b.coefs):
            s += a.norm[ia]*b.norm[ib]*ca*cb*\
                     Sxyz(a.exps[ia],a.shell,a.origin,
                     b.exps[ib],b.shell,b.origin)
    return s

