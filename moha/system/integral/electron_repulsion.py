import numpy as np
from moha.system.integral.nuclear_attraction import E,R
from ..auxiliary import gaussian_product_center

def gxyz(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,lmn4,D):
    ''' Evaluates kinetic energy integral between two Gaussians
         Returns a float.
         a,b,c,d:   orbital exponent on Gaussian 'a','b','c','d'
         lmn1,lmn2
         lmn3,lmn4: int tuple containing orbital angular momentum
                    for Gaussian 'a','b','c','d', respectively
         A,B,C,D:   list containing origin of Gaussian 'a','b','c','d'
    '''
    l1,m1,n1 = lmn1
    l2,m2,n2 = lmn2
    l3,m3,n3 = lmn3
    l4,m4,n4 = lmn4
    p = a+b # composite exponent for P (from Gaussians 'a' and 'b')
    q = c+d # composite exponent for Q (from Gaussians 'c' and 'd')
    alpha = p*q/(p+q)
    P = gaussian_product_center(a,A,b,B) # A and B composite center
    Q = gaussian_product_center(c,C,d,D) # C and D composite center
    RPQ = np.linalg.norm(P-Q)

    val = 0.0
    for t in range(l1+l2+1):
        for u in range(m1+m2+1):
            for v in range(n1+n2+1):
                for tau in range(l3+l4+1):
                    for nu in range(m3+m4+1):
                        for phi in range(n3+n4+1):
                            val += E(l1,l2,t,A[0],B[0],a,b) * \
                                   E(m1,m2,u,A[1],B[1],a,b) * \
                                   E(n1,n2,v,A[2],B[2],a,b) * \
                                   E(l3,l4,tau,C[0],D[0],c,d) * \
                                   E(m3,m4,nu ,C[1],D[1],c,d) * \
                                   E(n3,n4,phi,C[2],D[2],c,d) * \
                                   np.power(-1,tau+nu+phi) * \
                                   R(t+tau,u+nu,v+phi,0,\
                                       alpha,P[0]-Q[0],P[1]-Q[1],P[2]-Q[2],RPQ)

    val *= 2*np.power(np.pi,2.5)/(p*q*np.sqrt(p+q))
    return val

def electron_repulsion(a,b,c,d):
    '''Evaluates overlap between two contracted Gaussians
        Returns float.
        Arguments:
        a: contracted Gaussian 'a', BasisFunction object
        b: contracted Gaussian 'b', BasisFunction object
        c: contracted Gaussian 'b', BasisFunction object
        d: contracted Gaussian 'b', BasisFunction object
    '''
    eri = 0.0
    for ja, ca in enumerate(a.coefs):
        for jb, cb in enumerate(b.coefs):
            for jc, cc in enumerate(c.coefs):
                for jd, cd in enumerate(d.coefs):
                    eri += a.norm[ja]*b.norm[jb]*c.norm[jc]*d.norm[jd]*\
                             ca*cb*cc*cd*\
                             gxyz(a.exps[ja],a.shell,a.origin,\
                                b.exps[jb],b.shell,b.origin,\
                                c.exps[jc],c.shell,c.origin,\
                                d.exps[jd],d.shell,d.origin)
    return eri
