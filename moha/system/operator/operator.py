from moha.system.molecule import Molecule
from moha.system.basis.basis import Orbital
from moha.system.basis.basis_set import BasisSet

import numpy as np

__all__ = []

class Operator(np.ndarray):
    """General operator class.

    Methods
    -------
    __new__(cls)
        Build new instance.

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
        if len(args)==0:
            obj = super(Operator,cls).__new__(cls,0)
        elif len(args)==1:
            obj = np.asarray(args[0]).view(cls)
            if (len(set(obj.shape))!=1) or (1 not in set(obj.shape)):
                raise ValueError("Shape of the tensor must be general square.")
        else:
            raise ValueError("Number of arguments must be 0 or 1.")
        return obj

    def basis_transform(self,matrix):
        """Rotate orbitals with a transformation matrix.

        Parameters
        ----------
        matrix : np.ndarray(K, L)
            Transformation matrix.

        Raises
        ------
        TypeError
            If matrix is not a numpy array.
        ValueError
            If shape of matrix does not match up with the shape of the integrals.

        """
        if not isinstance(matrix, np.ndarray):
            raise TypeError("Transformation matrix must be given as a numpy array.")
        if matrix.ndim != 2:
            raise ValueError("Transformation matrix must be two-dimensional.")
        if matrix.shape[0] != self.shape[1]:
            raise ValueError(
                "Shape of the transformation matrix must match with the shape of the integrals."
            )
        integral = np.einsum("ij,ia->aj", self,matrix)
        integral = np.einsum("aj,jb->ab", integral,matrix)
        integral = integral.view(self.__class__)
        return integral
