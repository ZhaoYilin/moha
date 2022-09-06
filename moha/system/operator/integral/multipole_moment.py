import numpy as np

def MM(i,j,e,A,B,C,a,b):
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
    PC = P-C
    AB = A-B
    
    if i == j == e == 0:
        # base case
        return np.exp(-q*AB*AB) # K_AB

    elif i<0 or j<0 or e<0:
        return 0.

    elif i == 0 and j == 0 and e>0:
        # decrement index e
        return  PC*MM(i,j,e-1,A,B,C,a,b) +\
                1./(2*p)*i*MM(i-1,j,e,A,B,C,a,b) +\
                1./(2*p)*j*MM(i,j-1,e,A,B,C,a,b) +\
                1./(2*p)*(e-1)*MM(i,j,e-2,A,B,C,a,b)

    elif i == 0 and j>0 and e>0:
        # decrement index j
        return  PB*MM(i,j-1,e,A,B,C,a,b) +\
                1./(2*p)*i*MM(i-1,j-1,e,A,B,C,a,b) +\
                1./(2*p)*(j-1)*MM(i,j-2,e,A,B,C,a,b) +\
                1./(2*p)*e*MM(i,j-1,e-1,A,B,C,a,b)

    else:
        # decrement index i
        return  PA*MM(i-1,j,e,A,B,C,a,b) +\
                1./(2*p)*(i-1)*MM(i-2,j,e,A,B,C,a,b) +\
                1./(2*p)*j*MM(i-1,j-1,e,A,B,C,a,b) +\
                1./(2*p)*e*MM(i-1,j,e-1,A,B,C,a,b)
               

def MMxyz(a,lmn1,A,b,lmn2,B,lmn3,C):
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
    e,f,g = lmn3 # angular momentum
    MMx = MM(i,j,e,A[0],B[0],C[0],a,b) # X
    MMy = MM(k,l,f,A[1],B[1],C[1],a,b) # Y
    MMz = MM(m,n,g,A[2],B[2],C[2],a,b) # Z
    return np.power(np.pi/(a+b),0.5)*MMx*MMy*MMz
