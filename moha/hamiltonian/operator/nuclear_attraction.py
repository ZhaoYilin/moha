import itertools
import numpy as np

from moha.hamiltonian.operator.base import Operator

class NuclearAttraction(Operator):
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
        from moha.hamiltonian.integral.nuclear_attraction import NuclearAttraction
        obj = super().__new__(cls,*args)
        obj.integral = NuclearAttraction().contracted
        return obj

    def __call__(self, cga, cgb, mol):
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
        return self.integral(cga,cgb,mol)

    @classmethod
    def build(cls,molecule,basis_set):
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
        operator: Kinetic
            Return an operator instance.
        """
        from moha.molecule import Molecule
        from moha.basis.basis_set import BasisSet
        if not isinstance(molecule, Molecule):
            raise TypeError("molecule parameter must be a Molecule instance")
        if not isinstance(basis_set, BasisSet):
            raise TypeError("basis_set parameter must be a BasisSet instance")
        norb = len(basis_set)
        operator = cls(np.zeros((norb,norb)))
        for i, j in itertools.product(range(norb), repeat=2):
            if i <= j:
                operator[i,j] = operator.integral(basis_set[i],basis_set[j],molecule)
        operator = operator + operator.T - np.diag(operator.diagonal())
        return operator
