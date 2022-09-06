import numpy as np
from moha.system.operator.integral.overlap import S

def T(i,j,A,B,a,b):
    ''' Recursive definition of Hermite Gaussian coefficients.
        Returns a float.
        a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
        i,j: orbital angular momentum number on Gaussian 'a' and 'b'
        t: number nodes in Hermite (depends on type of integral,
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
        return (a-2*a**2*(PA**2+1./(2*p)))*S(0,0,A,B,a,b)

    elif i<0 or j<0 :
        return 0.

    elif i == 0 and j>0:
        # decrement index j
        return  PB*T(i,j-1,A,B,a,b) +\
                1./(2*p)*i*T(i-1,j-1,A,B,a,b) +\
                1./(2*p)*(j-1)*T(i,j-2,A,B,a,b) +\
                2*a*b/p*S(i,j,A,B,a,b) -\
                a/p*(j-1)*S(i,j-2,A,B,a,b) 

    else:
        # decrement index i
        return  PA*T(i-1,j,A,B,a,b) +\
                1./(2*p)*(i-1)*T(i-2,j,A,B,a,b) +\
                1./(2*p)*j*T(i-1,j-1,A,B,a,b) +\
                2*a*b/p*S(i,j,A,B,a,b) -\
                b/p*(i-1)*S(i-2,j,A,B,a,b) 


def Txyz(a,lmn1,A,b,lmn2,B):
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
    term0 = T(i,j,A[0],B[0],a,b)*S(k,l,A[1],B[1],a,b)*S(m,n,A[2],B[2],a,b)
    term1 = S(i,j,A[0],B[0],a,b)*T(k,l,A[1],B[1],a,b)*S(m,n,A[2],B[2],a,b)
    term2 = S(i,j,A[0],B[0],a,b)*S(k,l,A[1],B[1],a,b)*T(m,n,A[2],B[2],a,b)
    return np.power(np.pi/(a+b),1.5)*(term0+term1+term2)
