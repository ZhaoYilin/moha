from math import cos
from math import gcd
from math import pi
from math import sin
import sys

import numpy as np

__all__ = [
    "SymmetryOperator", "Identity", "Inversion", "Reflection", "Rotation",
    "ImproperRotation"
]


class SymmetryOperator:
    """ Base class of point group symmetry operator.
    """
    def __init__(self, symmetry_element, order):
        """ 
        """
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return np.allclose(self.matrix, other.matrix, rtol=1e-3)
        else:
            return False

    def __hash__(self):
        return hash(self.label)

    def __repr__(self):
        return self.label

    def __call__(self, mol):
        from moha.molecule import Molecule
        geometry = np.array(mol.geometry)
        geometry[:, 1:] = np.dot(self.matrix, geometry[:, 1:].T).T
        geometry = geometry.tolist()
        return Molecule.build(geometry, mol.pg)

    @property
    def matrix(self):
        """Returns the matrix representation of operator.

        Returns
        -------
        matrix : ndarray
            The matrix representation of operator.
        """
        return NotImplemented

    @property
    def homogeneous_matrix(self):
        """Returns the matrix representation of operator.

        Returns
        -------
        matrix : ndarray
            The matrix representation of operator.
        """
        if isinstance(self.matrix, type(NotImplemented)):
            return NotImplemented
        else:
            matrix = self.matrix
            homogeneous_matrix = np.zeros((4, 4))
            homogeneous_matrix[3, 3] = 1
            homogeneous_matrix[0:3, 0:3] = matrix
            return homogeneous_matrix


class Identity(SymmetryOperator):
    """ The identity operation.

    Attributes
    ----------
    order : int
        The smallest positive integer n, such that power n of the symmetry operator 
        equal to identity.

    label : str
        The symmetry operation label.

    Properties
    ----------
    matrix(self)
        Returns the matrix representation of operator.
    
    homogeneous_matrix(self)
        Returns the homogeneous matrix representation of operator.
    """
    def __init__(self, symmetry_element=1):
        """ Initialize the instance.

        Parameters
        ----------
        symmetry_element : int
            Identity element. 
        """
        self.order = 1
        self.symmetry_element = symmetry_element
        self.label = 'E'

    @property
    def matrix(self):
        """Returns the matrix representation of operator.

        Returns
        -------
        matrix : ndarray
            The matrix representation of operator.
        """
        matrix = np.identity(3)
        return matrix


class Inversion(SymmetryOperator):
    """ The Inversion symmetry operator.

    Attributes
    ----------
    order : int
        The smallest positive integer n, such that power n of the symmetry operator 
        equal to identity.

    label : str
        The symmetry operation label.
    """
    def __init__(self, symmetry_element=[0, 0, 0]):
        """ Initialize the instance.

        Parameters
        ----------
        symmetry_element : List[float,float,float]
            Coordinate of inversion center.
        """
        self.order = 2
        self.symmetry_element = symmetry_element
        self.label = 'i'

    def __call__(self, mol):
        """_summary_

        Parameters
        ----------
        mol : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """        
        from moha.molecule import Molecule
        natom = len(mol)
        geometry = np.array(mol.geometry)
        # Append ones column at the end
        homo_geometry = np.concatenate((geometry, np.ones((natom, 1))), axis=1)
        homo_geometry[:, 1:] = np.dot(self.homogeneous_matrix,
                                      homo_geometry[:, 1:].T).T
        # Remove ones column at the end
        geometry = homo_geometry[:, :-1]
        geometry = geometry.tolist()
        return Molecule.build(geometry, mol.pg)

    @property
    def matrix(self):
        """Returns the matrix representation of operator.

        Returns
        -------
        matrix : ndarray
            The matrix representation of operator.
        """
        if self.symmetry_element == [0, 0, 0]:
            matrix = -np.identity(3)
            return matrix
        else:
            raise ValueError("Center must be origin point.")

    @property
    def homogeneous_matrix(self):
        """Returns the matrix representation of operator.

        Returns
        -------
        matrix : ndarray
            The matrix representation of operator.
        """
        x, y, z = self.symmetry_element
        # Translation Matrix
        d = np.array([[1, 0, 0, x], [0, 1, 0, y], [0, 0, 1, z], [0, 0, 0, 1]])
        # Inverse of translation Matrix
        d_ = np.array([[1, 0, 0, -x], [0, 1, 0, -y], [0, 0, 1, -z],
                       [0, 0, 0, 1]])
        homogeneous_matrix = -np.identity(4)
        homogeneous_matrix[3, 3] = 1
        homogeneous_matrix = np.dot(d, np.dot(homogeneous_matrix, d_))
        return homogeneous_matrix


class Reflection(SymmetryOperator):
    """ The Reflection symmetry operator.

    This class can only create reflection planes that go through
    the origin.

    Attributes
    ----------
    order : int
        The smallest positive integer n, such that power n of the symmetry operator 
        equal to identity.

    label : str
        Label within the set {sigma_h, sigma_v, sigma_d}
        h : horizontal; v : vertical; d : dihedral

    symmetry_element : List[float, float, float]  
        Normal vector of reflection plane.
    """
    def __init__(self, symmetry_element):
        """ Initialize the instance.

        Parameters
        ----------
        symmetry_element : List[float,float,float]
            Normal vector of reflection plane.
        """
        self.order = 2
        self.label = 'sigma'
        self.symmetry_element = symmetry_element

    @property
    def matrix(self):
        """Returns the matrix representation of operator.

        Returns
        -------
        matrix : ndarray
            The matrix representation of operator.
        """
        norm = self.symmetry_element / np.linalg.norm(self.symmetry_element)
        matrix = np.eye(3) - 2 * np.outer(norm, norm)
        return matrix


class Rotation(SymmetryOperator):
    """Rotation symmetry operator.

    This class can only create axis of rotation that go through
    the origin.

    Attributes
    ----------
    order : int
        The smallest positive integer n, such that power n of the symmetry operator 
        equal to identity.

    n : int
        The rotation order of fold.

    axis : List[float, float, float]
        The axis of rotation.
    """
    def __init__(self, symmetry_element, order, k=1):
        """ Initialize the instance.

        Parameters
        ----------
        symmetry_element : List[float,float,float]
            Vector of the proper rotation axis. 

        n : int
            The rotation order of fold.

        k : int
            The rotation order of fold.
        """
        self.symmetry_element = symmetry_element
        self.order = order
        self.k = k
        self.label = "C" + "_{" + str(order) + "}" + "^{" + str(k) + "}"

    def __repr__(self):
        """_summary_

        Returns
        -------
        _type_
            _description_
        """        
        if self.order == sys.maxsize:
            return "C_{oo}"
        else:
            return self.label

    @property
    def matrix(self):
        """Rotates a point around the axis of rotation for a angle of
        2 * pi / self.fold

        Parameters
        ----------
        coordinate : Tuple[float, float, float]

        Returns
        -------
        : Tuple[float, float, float]

        """
        theta = 2 * pi * self.k / self.order
        norm = self.symmetry_element / np.linalg.norm(self.symmetry_element)
        x, y, z = norm
        xx = x**2 * (1 - cos(theta)) + cos(theta)
        xy = x * y * (1 - cos(theta)) - z * sin(theta)
        xz = x * z * (1 - cos(theta)) + y * sin(theta)
        yx = x * y * (1 - cos(theta)) + z * sin(theta)
        yy = y**2 * (1 - cos(theta)) + cos(theta)
        yz = y * z * (1 - cos(theta)) - x * sin(theta)
        zx = x * z * (1 - cos(theta)) - y * sin(theta)
        zy = z * y * (1 - cos(theta)) + x * sin(theta)
        zz = z**2 * (1 - cos(theta)) + cos(theta)
        matrix = np.array([[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]])
        return matrix

    @property
    def expansion(self):
        """_summary_

        Returns
        -------
        _type_
            _description_
        """ 
        C_list = []
        for k in range(self.order)[1:]:
            C = self.__class__(self.symmetry_element, self.order, k)
            C_list.append(C)
        return C_list


class ImproperRotation(SymmetryOperator):
    """Improper rotation symmetry operator.

    This class can only create axis of improper rotation that goes
    through the origin. As the reflection plane of the improper rotation
    will go through the origin, it is therefore equivalent to rotation
    and inversion. This operator will instead rotate and invert as the
    code is simpler and does not need the creation of the householder
    matrix used in reflection operations.

    Attributes
    ----------
    order : int
        The smallest positive integer n, such that power n of the symmetry operator 
        equal to identity.

    fold : float
        The fraction of 2 * pi angle a coordinate is rotated.

    vector : Tuple[float, float, float]
        The axis of improper rotation.
    """
    def __init__(self, symmetry_element, order, k=1):
        """ Initialize the instance.

        Parameters
        ----------
        symmetry_element : List[float,float,float]
            Vector of the proper rotation axis. 

        n : int
            The rotation order of fold.

        k : int
            The rotation order of fold.
        """
        self.symmetry_element = symmetry_element
        self.order = order
        self.k = k
        self.label = "S" + "_{" + str(order) + "}" + "^{" + str(k) + "}"

    def __repr__(self):
        if self.order == sys.maxsize or self.order == 2 * sys.maxsize:
            return "S_{oo}"
        else:
            return self.label

    @property
    def matrix(self):
        """Rotates a point around the axis of rotation for a angle of
        2 * pi / self.fold and then inverts.

        Parameters
        ----------
        coordinate : Tuple[float, float, float]

        Returns
        -------
        : Tuple[float, float, float]

        """
        R_matrix = Rotation(self.symmetry_element, self.order, self.k).matrix
        sigma_matrix = Reflection(self.symmetry_element).matrix
        matrix = np.dot(R_matrix, sigma_matrix)
        return matrix
