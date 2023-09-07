import itertools
import numpy as np
from multiprocessing import Pool

from moha.hamiltonian.operator.base import Operator

class ElectronRepulsion(Operator):
    """Electron repulsion operator class.
    
    Methods
    -------
    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __new__(cls,*args):
        """Generate new instance.
        Parameters
        ----------
        *args : null or iterator
            Null for empty instance.
            Iterator for assign np.array.
        """
        from moha.hamiltonian.integral.electron_repulsion import ElectronRepulsion
        obj = super().__new__(cls,*args)
        obj.integral = ElectronRepulsion().contracted
        return obj

    def __call__(self, cga, cgb, cgc, cgd):
        """Evaluates overlap between two orbitals.

        Parameters
        ----------
        orb1: Orbital
            The first orbital.

        orb2: Orbital
            The second orbital.

        orb3: Orbital
            The third orbital.

        orb4: Orbital
            The fourth orbital.

        Return
        ------
        vlaue: float
            Integral value.
        """
        return self.integral(cga,cgb,cgc,cgd)

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
        operator: Operator
            Return an operator instance.
        """
        from moha.basis.basis_set import BasisSet
        from moha.hamiltonian.auxiliary import eint
        if not isinstance(basis_set, BasisSet):
            raise TypeError("basis_set parameter must be a BasisSet instance")

        norb = len(basis_set)
        integral_list = []    
        operator = cls(np.zeros((norb,norb,norb,norb)))

        # Calculate integral list with 8-fold permutational symmetry
        indices = []
        bases = []
        n = 0
        for a, b, c, d in itertools.product(range(norb), repeat=4):
            if not (a < b or c < d or a < c or (a == c and b < d)):
                indices.append((a,b,c,d))
                bases.append((basis_set[a],basis_set[b],basis_set[c],basis_set[d]))
        pool = Pool(operator.ncpu)
        integral_list = pool.starmap(operator.integral, bases)
        pool.close()

        # Build integral tensor from list 
        for a, b, c, d in itertools.product(range(norb), repeat=4):
            operator[a,b,c,d] = integral_list[eint(a,b,c,d)]

        return operator

    @property
    def coulomb(self):
        """Coulomb integral over spatial orbitals.
        Jij = (ii|jj)

        Returns
        -------
        Jij : np.ndarray
            Coulomb integral.        
        """
        Jij = np.einsum("iijj->ij", self)
        return Jij 

    @property
    def exchange(self):
        """Exchange integral over spatial orbitals.        
        Kij = (ij|ij)

        Returns
        -------
        Kij : np.ndarray
            Exchange integral.        
        """
        Kij = np.einsum("ijij->ij", self)
        return Kij

    def basis_transform(self,C):
        """Transform Eri from atomic orbitals to molecular orbitals.

        Parameters
        ----------
        C : np.ndarray
            Coefficient matrix.

        Raises
        ------
        TypeError
            If coefficient matrix is not a numpy array.
        ValueError
            If shape of coefficient matrix does not match up with the shape of
            the integrals.
        """
        if not isinstance(C, np.ndarray):
            raise TypeError("Coefficient matrix must be given as a numpy array.")
        if C.ndim != 2:
            raise ValueError("Coefficient matrix must be two-dimensional.")
        if C.shape[0] != self.shape[1]:
            raise ValueError(
                "Shape of the coefficient matrix must match with the shape of Eri."
            )
        Eri = np.einsum("ijkl,ia->ajkl", self, C)
        Eri = np.einsum("ajkl,jb->abkl", Eri, C)
        Eri = np.einsum("abkl,kc->abcl", Eri, C)
        Eri = np.einsum("abcl,ld->abcd", Eri, C)
        Eri = Eri.view(self.__class__)
        return Eri
