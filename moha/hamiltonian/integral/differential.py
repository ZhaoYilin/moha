import numpy as np
import itertools

class Differential:
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

    def contracted(self, cga, cgb, order):
        """Evaluate integral over two contracted GTOs.

        Parameters
        ----------
        cga: ContractedGaussian
            The first contracted gaussian orbital.

        order : List[int, int, int]
            Differentail order of x y and z direction.

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
            result += contraction*self.primitive(pga, pgb, order)
        return result

    def primitive(self, pga, pgb, order):
        """Evaluate integral over two primitive gaussian orbitals.

        Parameters
        ----------
        pga: PrimitiveGaussian
            The first primitive gaussian orbital.

        order : List[int, int, int]
            Differentail order of x y and z direction.

        pgb: PrimitiveGaussian
            The second primitive gaussian orbital.
        
        Return
        ------
        result : float
            Integral value.
        """
        result = 1
        for r in range(3):
            a = pga.exponent
            b = pgb.exponent
            A = pga.origin[r]
            B = pgb.origin[r]
            i = pga.shell[r]
            j = pgb.shell[r]
            e = order[r]
            result *= self.D1d(A, i, a, B, j, b, e)
        return result

    def D1d(self, A, i, a, B, j, b, e):
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
        if i<0 or j<0 or e<0:
            return 0.0
        elif e == 1:
            return 2*a*self.overlap.S1d(A,i+1,a,B,j,b)-i*self.overlap.S1d(A,i-1,a,B,j,b)
        elif e == 2:
            return 4*a**2*self.overlap.S1d(A,i+2,a,B,j,b) -\
                2*a*(2*i+1)*self.overlap.S1d(A,i,a,B,b) +\
                i*(i-1)*self.overlap.S1d(A,i-2,a,B,j,b)
        else:
            return 2*a*self.D1d(A,i+1,a,B,j,b,e-1) - i*self.D1d(A,i-1,a,B,j,b,e-1)

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

    D = Differential()
    d17 = D.contracted(cga,cgb,[1,0,0])
    print(d17)
    #print(np.isclose(t17,-0.167203))
