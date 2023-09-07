import os
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
        else:
            raise ValueError("Number of arguments must be 0 or 1.")

        obj.ncpu = os.cpu_count()
        return obj

    def basis_transform(self,matrix):
        """Rotate orbitals with a transformation matrix.

        Parameters
        ----------
        matrix : np.ndarray(K, L)
            Transformation matrix.

        """
        integral = np.einsum("ij,ia->aj", self,matrix)
        integral = np.einsum("aj,jb->ab", integral,matrix)
        integral = integral.view(self.__class__)
        return integral
