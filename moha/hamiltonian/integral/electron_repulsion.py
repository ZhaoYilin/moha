import itertools
import numpy as np

from math import dist, exp, pi, sqrt
from moha.hamiltonian.integral.boys import *
from moha.basis.gauss import PrimitiveGaussian

class ElectronRepulsion:
    """The Obara-Saika scheme for three-dimensional nuclear attraction integral over
    primitive Gaussian orbitals.

    Attributes
    ----------

    Methods
    -------
    __init__(self)
        Initialize the instance.

    contracted(self, cga, cgb, cgc, cgd)
        Evaluate integral over contracted GTOs.

    target(self, N, pga, pgb, pgc, pgd)
        Obtain the target integral in the Obara-Saika recurrence relations.

    recurrence(self, r, N, pga, pgb, pgc, pgd, pga_1, pga_2, pgb_1, pgc_1, pgd_1)
        Run the recurrence.

    new_gaussian(self, r, pga, pgb, pgc, pgd)
        Generate all the new primitive GTOs in the Obara-Saika recurrence relations.
    """
    def __init__(self):
        """Initialize the instance.
        """
        self.alpha = 0
        self.R = ()
        self.boys_dict = {}

    def contracted(self, cga, cgb, cgc, cgd):
        """Evaluate integral over contracted GTOs.

        Parameters
        ----------
        cga: ContractedGaussian
            The first contracted gaussian orbital.

        cgb: ContractedGaussian
            The second contracted gaussian orbital.

        cgc: ContractedGaussian
            The third contracted gaussian orbital.

        cgd: ContractedGaussian
            The fourth contracted gaussian orbital.

        Return
        ------
        result : float
            Integral value.
        """
        l_total = sum(cga.shell) + sum(cgb.shell) + sum(cgc.shell) + sum(cgd.shell)

        A = np.array(cga.origin)
        B = np.array(cgb.origin)
        C = np.array(cgc.origin)
        D = np.array(cgd.origin)
        RAB = dist(A,B)
        RCD = dist(C,D)

        primitives_a = cga.expansion
        primitives_b = cgb.expansion
        primitives_c = cgc.expansion
        primitives_d = cgd.expansion

        result = 0.0
        for pga, pgb, pgc, pgd in itertools.product(primitives_a, primitives_b, primitives_c, primitives_d):
            ca = pga.coefficient
            cb = pgb.coefficient
            cc = pgc.coefficient
            cd = pgd.coefficient
            na = pga.norm
            nb = pgb.norm
            nc = pgc.norm
            nd = pgd.norm
            contraction = ca * cb * cc * cd * na * nb * nc * nd

            a = pga.exponent
            b = pgb.exponent
            c = pgc.exponent
            d = pgd.exponent
            p = a + b
            q = c + d
            mu = (a*b)/(a+b)
            nu = (c*d)/(c+d)
            self.alpha = (p * q) / (p + q)

            P = (a*A+b*B)/(a+b)
            Q = (c*C+d*D)/(c+d)
            self.R = (p*P+q*Q)/(p+q)
            RPQ = dist(P,Q)

            # Build boys function F_{N}(x)
            Kab = exp(-mu*RAB**2)
            Kcd = exp(-nu*RCD**2)
            boys_pre_factor = (2*pi**(5/2))/(p*q*sqrt(p+q))*Kab*Kcd
            N = l_total
            x = self.alpha*RPQ**2
            boys_function = boys(l_total, x)
            Theta_N_0000_0000_0000 = boys_pre_factor * boys_function
            self.boys_dict = {l_total: Theta_N_0000_0000_0000}

            while N >= 1:
                boys_function = boys_recursion(N, x, boys_function)
                N -= 1
                Theta_N_0000_0000_0000 = boys_pre_factor * boys_function
                self.boys_dict[N] = Theta_N_0000_0000_0000
            result += contraction * self.target(0, pga, pgb, pgc, pgd)
        return result

    def target(self, N, pga, pgb, pgc, pgd):
        """Obtain the target integral in the Obara-Saika recurrence relations.

        Parameters
        ----------
        N : int
            Order of the boys function F_{N}(x).

        pga : PrimitiveGaussian
            The primitive gaussian orbital.

        pgb : PrimitiveGaussian
            The primitive gaussian orbital.

        pgc : PrimitiveGaussian
            The primitive gaussian orbital.

        pgd : PrimitiveGaussian
            The primitive gaussian orbital.
        """
        if pga.shell[0] > 0:
            return self.recurrence(0, N, pga, pgb, pgc, pgd, *self.new_gaussian(0, pga, pgb, pgc, pgd))
        elif pga.shell[1] > 0:
            return self.recurrence(1, N, pga, pgb, pgc, pgd, *self.new_gaussian(1, pga, pgb, pgc, pgd))
        elif pga.shell[2] > 0:
            return self.recurrence(2, N, pga, pgb, pgc, pgd, *self.new_gaussian(2, pga, pgb, pgc, pgd))
        elif pgb.shell[0] > 0:
            return self.recurrence(0, N, pgb, pga, pgd, pgc, *self.new_gaussian(0, pgb, pga, pgd, pgc))
        elif pgb.shell[1] > 0:
            return self.recurrence(1, N, pgb, pga, pgd, pgc, *self.new_gaussian(1, pgb, pga, pgd, pgc))
        elif pgb.shell[2] > 0:
            return self.recurrence(2, N, pgb, pga, pgd, pgc, *self.new_gaussian(2, pgb, pga, pgd, pgc))
        elif pgc.shell[0] > 0:
            return self.recurrence(0, N, pgc, pgd, pga, pgb, *self.new_gaussian(0, pgc, pgd, pga, pgb))
        elif pgc.shell[1] > 0:
            return self.recurrence(1, N, pgc, pgd, pga, pgb, *self.new_gaussian(1, pgc, pgd, pga, pgb))
        elif pgc.shell[2] > 0:
            return self.recurrence(2, N, pgc, pgd, pga, pgb, *self.new_gaussian(2, pgc, pgd, pga, pgb))
        elif pgd.shell[0] > 0:
            return self.recurrence(0, N, pgd, pgc, pgb, pga, *self.new_gaussian(0, pgd, pgc, pgb, pga))
        elif pgd.shell[1] > 0:
            return self.recurrence(1, N, pgd, pgc, pgb, pga, *self.new_gaussian(1, pgd, pgc, pgb, pga))
        elif pgd.shell[2] > 0:
            return self.recurrence(2, N, pgd, pgc, pgb, pga, *self.new_gaussian(2, pgd, pgc, pgb, pga))
        else:
            return self.boys_dict[N]

    def recurrence(self, r, N, pga, pgb, pgc, pgd, pga_1, pga_2, pgb_1, pgc_1, pgd_1):
        """Run the recurrence.

        Parameters
        ----------
        r : int
            Cartesian index 0, 1, 2. 

        N : int
            Order of the boys function F_{N}(x).

        pga : PrimitiveGaussian
            The first primitive gaussian orbital with angular momentum i.

        pgb : PrimitiveGaussian
            The second primitive gaussian orbital with angular momentum j.

        pgc : PrimitiveGaussian
            The third primitive gaussian orbital with angular momentum k.

        pgd : PrimitiveGaussian
            The fourth primitive gaussian orbital with angular momentum l.

        pga_1 : PrimitiveGaussian
            The first primitive gaussian orbital with angular momentum i-1.

        pga_2 : PrimitiveGaussian
            The first primitive gaussian orbital with angular momentum i-2.

        pgb_1 : PrimitiveGaussian
            The second primitive gaussian orbital with angular momentum j-1.

        pgc_1 : PrimitiveGaussian
            The third primitive gaussian orbital with angular momentum k-1.

        pgd_1 : PrimitiveGaussian
            The fourth primitive gaussian orbital with angular momentum l-1.

        Return
        ------
        result : float
            Integral value.
        """
        term1 = term2 = term3 = term4 = term5 = term6 = term7 = term8 = 0

        a = pga.exponent
        b = pgb.exponent
        c = pgc.exponent
        d = pgd.exponent
        p = a + b
        q = c + d

        A = np.array(pga.origin)
        B = np.array(pgb.origin)
        C = np.array(pgc.origin)
        D = np.array(pgd.origin)
        P = (a*A+b*B)/(a+b)
        Q = (c*C+d*D)/(c+d)

        XPA = P - A
        XPQ = P - Q

        if XPA[r] != 0:
            term1 = XPA[r] * self.target(N, pga_1, pgb, pgc, pgd)
        if XPQ[r] != 0:
            term2 = self.alpha/p*XPQ[r] * self.target(N+1, pga_1, pgb, pgc, pgd)
        if pga_1.shell[r] > 0:
            term3 = pga_1.shell[r] * (1 / (2 * p)) * self.target(N, pga_2, pgb, pgc, pgd)
            term4 = pga_1.shell[r] * (self.alpha / (2 * p ** 2)) * self.target(N+1, pga_2, pgb, pgc, pgd)
        if pgb.shell[r] > 0:
            term5 = pgb.shell[r] * (1 / (2 * p)) * self.target(N, pga_1, pgb_1, pgc, pgd)
            term6 = pgb.shell[r] * (self.alpha / (2 * p ** 2)) * self.target(N+1, pga_1, pgb_1, pgc, pgd)
        if pgc.shell[r] > 0:
            term7 = pgc.shell[r] * (1 / (2 * (p + q))) * self.target(N+1, pga_1, pgb, pgc_1, pgd)
        if pgd.shell[r] > 0:
            term8 = pgd.shell[r] * (1 / (2 * (p + q))) * self.target(N+1, pga_1, pgb, pgc, pgd_1)

        return term1 - term2 + term3 - term4 + term5 - term6 + term7 + term8

    def new_gaussian(self, r, pga, pgb, pgc, pgd):
        """Generate all the new primitive GTOs in the Obara-Saika recurrence relations.

        Parameters
        ----------
        r : int
            Cartesian index 0, 1, 2. 

        N : int
            Order of the boys function F_{N}(x).

        pga : PrimitiveGaussian
            The first primitive gaussian orbital.

        pgb : PrimitiveGaussian
            The second primitive gaussian orbital.

        pgc : PrimitiveGaussian
            The third primitive gaussian orbital.

        pgd : PrimitiveGaussian
            The fourth primitive gaussian orbital.

        Return
        ------
        result : Tuple(pg,pg,pg,pg,pg)
            Tuple of 5 PrimitiveGaussian orbital instance. 
        """
        ca = pga.coefficient
        cb = pgb.coefficient
        cc = pgc.coefficient
        cd = pgd.coefficient

        a = pga.exponent
        b = pgb.exponent
        c = pgc.exponent
        d = pgd.exponent

        A = pga.origin
        B = pgb.origin
        C = pgc.origin
        D = pgd.origin

        ix,iy,iz = pga.shell
        jx,jy,jz = pgb.shell
        kx,ky,kz = pgc.shell
        lx,ly,lz = pgd.shell

        if r == 0:
            pga_1 = PrimitiveGaussian(ca, A, (ix - 1, iy, iz), a)
            pga_2 = PrimitiveGaussian(ca, A, (ix - 2, iy, iz), a)
            pgb_1 = PrimitiveGaussian(cb, B, (jx - 1, jy, jz), b)
            pgc_1 = PrimitiveGaussian(cc, C, (kx - 1, ky, kz), c)
            pgd_1 = PrimitiveGaussian(cd, D, (lx - 1, ly, lz), d)
            return pga_1, pga_2, pgb_1, pgc_1, pgd_1
        elif r == 1:
            pga_1 = PrimitiveGaussian(ca, A, (ix, iy-1, iz), a)
            pga_2 = PrimitiveGaussian(ca, A, (ix, iy-2, iz), a)
            pgb_1 = PrimitiveGaussian(cb, B, (jx, jy-1, jz), b)
            pgc_1 = PrimitiveGaussian(cc, C, (kx, ky-1, kz), c)
            pgd_1 = PrimitiveGaussian(cd, D, (lx, ly-1, lz), d)
            return pga_1, pga_2, pgb_1, pgc_1, pgd_1
        elif r == 2:
            pga_1 = PrimitiveGaussian(ca, A, (ix, iy, iz-1), a)
            pga_2 = PrimitiveGaussian(ca, A, (ix, iy, iz-2), a)
            pgb_1 = PrimitiveGaussian(cb, B, (jx, jy, jz-1), b)
            pgc_1 = PrimitiveGaussian(cc, C, (kx, ky, kz-1), c)
            pgd_1 = PrimitiveGaussian(cd, D, (lx, ly, lz-1), d)
            return pga_1, pga_2, pgb_1, pgc_1, pgd_1
