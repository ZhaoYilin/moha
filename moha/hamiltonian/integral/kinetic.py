import numpy as np
import itertools

class Kinetic:
    """The Obara-Saika scheme for kinetic energy integral.

    Methods
    -------
    __init__(self)
        Initialize the instance.

    contracted(self, cga, cgb)
        Evaluates integral over contracted GTOs.

    primitive(self, pga, pgb)
        Evaluates integral over primitive GTOs.

    target(self, r, pga, pgb)
        Obtain the target integral in the Obara-Saika recurrence relations.

    recurrence(self, r, pga, pgb, pga_1, pga_2, pgb_1)
        Run the recurrence.

    new_gaussian(self, r, pga, pgb):
        Generate all the new primitive GTOs in the Obara-Saika recurrence relations.
    """
    def __init__(self):
        """Initialize the instance.
        """
        from moha.hamiltonian.integral.overlap import Overlap
        self.overlap = Overlap()

    def contracted(self, cga, cgb):
        """Evaluate integral over two contracted GTOs.

        Parameters
        ----------
        cga: ContractedGaussian
            The first contracted gaussian orbital.

        cgb: ContractedGaussian
            The second contracted gaussian orbital.

        Return
        ------
        result : float
            Integral value.
        """
        primitives_a = cga.expansion
        primitives_b = cgb.expansion
        result = 0.0
        for pga, pgb in itertools.product(primitives_a, primitives_b):
            ca = pga.coefficient
            cb = pgb.coefficient
            na = pga.norm
            nb = pgb.norm
            contraction = ca * cb * na * nb
            result += contraction*self.primitive(pga,pgb)
        return result

    def primitive(self, pga, pgb):
        """Evaluate integral over two primitive gaussian orbitals.

        Parameters
        ----------
        pga: PrimitiveGaussian
            The first primitive gaussian orbital.

        pgb: PrimitiveGaussian
            The second primitive gaussian orbital.
    
        Return
        ------
        result : float
            Integral value.
        """
        a = pga.exponent
        b = pgb.exponent
        Ax,Ay,Az = pga.origin
        Bx,By,Bz = pgb.origin
        i,k,m = pga.shell
        j,l,n = pgb.shell

        Sij = self.overlap.S1d(Ax,i,a,Bx,j,b)
        Skl = self.overlap.S1d(Ay,k,a,By,l,b)
        Smn = self.overlap.S1d(Az,m,a,Bz,n,b)


        Tij = self.T1d(Ax,i,a,Bx,j,b)
        Tkl = self.T1d(Ay,k,a,By,l,b)
        Tmn = self.T1d(Az,m,a,Bz,n,b)

        Tab = Tij*Skl*Smn+Sij*Tkl*Smn+Sij*Skl*Tmn
        return Tab

    def T1d(self, A, i, a, B, j, b):
        """Obtain the target integral in the Obara-Saika recurrence relations.

        Parameters
        ----------
        A : float
            Origin of the basis for one direction.

        i : int
            Angular momentum for one direction.

        a : float
            Primitive Gaussian exponent for one direction.

        B : float
            Origin of the basis for one direction.

        j : int
            Angular momentum for one direction.

        b : float
            Primitive Gaussian exponent for one direction.

        Return
        ------
        result : float
            Integral value for one direction.
        """
        p = a + b
        P = (1./p)*(a*A + b*B)

        if i<0 or j<0:
            return 0.
        elif i>0:
            # decrement index i
            return  (P-A)*self.T1d(A,i-1,a,B,j,b) +\
                    1./(2*p)*(i-1)*self.T1d(A,i-2,a,B,j,b) +\
                    1./(2*p)*j*self.T1d(A,i-1,a,B,j-1,b) +\
                    2*a*b/p*self.overlap.S1d(A,i,a,B,j,b) -\
                    b/p*(i-1)*self.overlap.S1d(A,i-2,a,B,j,b)
        elif j>0:
            # decrement index j
            return  (P-B)*self.T1d(A,i,a,B,j-1,b) +\
                    1./(2*p)*i*self.T1d(A,i-1,a,B,j-1,b) +\
                    1./(2*p)*(j-1)*self.T1d(A,i,a,B,j-2,b) +\
                    2*a*b/p*self.overlap.S1d(A,i,a,B,j,b) -\
                    a/p*(j-1)*self.overlap.S1d(A,i,a,B,j-2,b)
        else:
            # boundary conditon i==j==0
            T00 = (a-2*a**2*((P-A)**2+1./(2*p)))*self.overlap.S1d(A,0,a,B,0,b)
            return T00

if __name__ == '__main__':
    # Coordinate of H2O molecule
    H2O = [[0., 1.43233673, -0.96104039],
    [0., -1.43233673, -0.96104039],
    [0., 0., 0.24026010]]

    # Primitive contraction coefficients
    PrimCoeff = np.array([[0.1543289673, 0.5353281423, 0.4446345422],
    [0.1543289673, 0.5353281423, 0.4446345422],
    [0.1543289673, 0.5353281423, 0.4446345422],
    [-0.09996722919, 0.3995128261, 0.7001154689],
    [0.155916275, 0.6076837186, 0.3919573931],
    [0.155916275, 0.6076837186, 0.3919573931],
    [0.155916275, 0.6076837186, 0.3919573931]])

    # Orbital exponents
    OrbCoeff = np.array([[3.425250914, 0.6239137298, 0.168855404],
    [3.425250914, 0.6239137298, 0.168855404],
    [130.7093214, 23.80886605, 6.443608313],
    [5.033151319, 1.169596125, 0.38038896],
    [5.033151319, 1.169596125, 0.38038896],
    [5.033151319, 1.169596125, 0.38038896],
    [5.033151319, 1.169596125, 0.38038896]])

    # H1s, H2s, O1s, O2s, O2px , O2py, O2p
    FCenter = [H2O[0], H2O[1], H2O[2], H2O[2], H2O[2], H2O[2], H2O[2]]
    CartAng = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0],
    [1, 0, 0], [0, 1, 0], [0, 0, 1]]

    from moha.basis.gauss import ContractedGaussian
    cga = ContractedGaussian(PrimCoeff[0],FCenter[0],CartAng[0],OrbCoeff[0])
    cgb = ContractedGaussian(PrimCoeff[6],FCenter[6],CartAng[6],OrbCoeff[6])

    T = Kinetic()
    t17 = T.contracted(cga,cgb)
    print(np.isclose(t17,-0.167203))
