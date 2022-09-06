import numpy as np
from moha.system.operator.integral.overlap import S

def D(i,j,e,A,B,a,b):
    ''' Recursive definition of Hermite Gaussian coefficients.
        Returns a float.
        a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
        i,j: orbital angular momentum number on Gaussian 'a' and 'b'
        t: number nodes in Hermite (depends on type of integral, 
           e.g. always zero for overlap integrals)
        Qx: distance between origins of Gaussian 'a' and 'b'
    '''
    if e < 1:
        raise ValueError('e should be integer larger than zero')
    elif e == 1:
        return 2*a*S(i+1,j,A,B,a,b)-i*S(i-1,j,A,B,a,b)
    elif e == 2:
        return 4*a**2*S(i+2,j,A,B,a,b) -\
               2*a*(2*i+1)*S(i,j,A,B,a,b) +\
               i*(i-1)*S(i-2,j,A,B,a,b)
    else:
        return 2*a*D(i+1,j,e-1,A,B,a,b) - i*D(i-1,j,e-1,A,B,a,b)

def Dxyz(a,lmn1,A,b,lmn2,B,es):
    ''' Evaluates overlap integral between two Gaussians
        Returns a float.
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
    e,g,f = es # shell angular momentum on Gaussian 'b'
    Dx = D(i,j,e,A[0],B[0],a,b) # X
    Dy = D(k,l,g,A[1],B[1],a,b) # Y
    Dz = D(m,n,f,A[2],B[2],a,b) # Z
    return np.power(np.pi/(a+b),1.5)*Dx*Dy*Dz
