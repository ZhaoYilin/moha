import numpy as np

class U1(object):
    """Particle number symmetry.

    Attributes
    ----------
    n: int
        Number of electrons.

    Properties
    ----------
    parity: int
        Parity of the system.
    """
    def __init__(self, n=0):
        """Initialize the instance.

        Parameters
        ----------
        n: int
            Number of electrons.
        """
        self.n = n

    def __add__(self, other):
        return self.__class__(self.n + other.n)

    def __sub__(self, other):
        return self.__class__(self.n - other.n)

    def __neg__(self):
        return self.__class__(-self.n)

    def __eq__(self, other):
        return self.n == other.n

    def __lt__(self, other):
        return self.n < other.n

    def __repr__(self):
        return "< U1(%d) >" % (self.n)

    @property
    def parity(self):
        return self.n % 2

