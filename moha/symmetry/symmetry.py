import numpy as np

__all__ = ['Symmetry']

class Symmetry:
    """Globle symmetries of the system: U(1) particle numer, SU(2) total 
    electronic spin, P abelian point group symmetry.

    Attributes
    ----------
    n: int
        Particle number.

    ms2: int
        The spin polarization of the system

    ipg: int
        The irreducible representation of the point group symmetry.

    Methods
    -------
    __init__(self, n=0, ms2=0, ipg=0)
        Initialize the instance.
    """

    def __init__(self, n=0, ms2=0, ipg=0):
        """Initialize the instance.

        Parameters
        ----------
        n: int
            Number of electrons

        ms2: int
            The spin polarization of the system

        ipg: int
            The irreducible representation of the point group symmetry.
        """
        self.n = n
        self.ms2 = ms2
        self.ipg = ipg

    def __add__(self, other):
        return self.__class__(self.n + other.n, self.ms2 + other.ms2, self.ipg ^ other.ipg)

    def __sub__(self, other):
        return self.__class__(self.n - other.n, self.ms2 - other.ms2, self.ipg ^ other.ipg)

    def __neg__(self):
        return self.__class__(-self.n, -self.ms2, self.ipg)

    def __eq__(self, other):
        return self.n == other.n and self.ms2 == other.ms2 and self.ipg == other.ipg

    def __lt__(self, other):
        return (self.n, self.ms2, self.ipg) < (other.n, other.ms2, other.ipg)

    def __repr__(self):
        if self.ms2 % 2 == 1:
            return "< N={0:d} SZ={1:d}/2 PG={2:d} >".format(self.n, self.ms2, self.ipg)
        else:
            return "< N={0:d} SZ={1:d} PG={2:d} >".format(self.n, self.ms2//2, self.ipg)
