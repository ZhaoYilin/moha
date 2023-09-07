import itertools
import numpy as np

from moha.hamiltonian.operator.base import Operator

class Kinetic(Operator):
    """Overlap operator class.
    
    Methods
    -------
    """
    def __new__(cls,*args):
        """Generate new instance.
        Parameters
        ----------
        *args : null or iterator
            Null for empty instance.
            Iterator for assign np.array.
        """
        from moha.hamiltonian.integral.kinetic import Kinetic
        obj = super().__new__(cls,*args)
        obj.integral = Kinetic().contracted
        return obj

    def __call__(self, cga, cgb):
        """Evaluates overlap between two orbitals.

        Parameters
        ----------
        orb1: Orbital
            The first orbital.

        orb2: Orbital
            The second orbital.

        Return
        ------
        vlaue: float
            Integral value.
        """
        return self.integral(cga,cgb)

    @classmethod
    def build(cls,basis_set):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : BasisSet
            Hatree Fock basis set instance.

        Raises
        ------
        TypeError
            If basis_set parameter is not a BasisSet instance.

        Return
        ------
        operator: Overlap
            Return an operator instance.
        """
        from moha.basis.basis_set import BasisSet
        if not isinstance(basis_set, BasisSet):
            raise TypeError("basis_set parameter must be a BasisSet instance.")
        norb = len(basis_set)    
        operator = cls(np.zeros((norb,norb)))
        for i, j in itertools.product(range(norb), repeat=2):
            if i <= j:
                operator[i,j] = operator.integral(basis_set[i],basis_set[j])
        operator = operator + operator.T - np.diag(operator.diagonal())
        return operator
