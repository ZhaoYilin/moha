from enum import Enum, auto
from abc import abstractmethod,abstractproperty
import numpy as np

class OperatorNames(Enum):
    """Operator Names.

    Members
    ----------
    I : enum.auto
        Identity operator.
    
    S : enum.auto
        Overlap operator.
    
    T : enum.auto
        Kinetic energy operator.
    
    V : enum.auto
        Nuclear attraction operator.
    
    Hcore : enum.auto
        Core Hamiltonian operator.
    
    MM : enum.auto
        Multipole moment operator.
    
    Diff : enum.auto
        Differential operator.
    
    LM : enum.auto
        Linear momentum operator.
    
    AM : enum.auto
        Angular momentum operator.
    
    Eri : enum.auto
        Electron repulsion operator.
    """
    I = auto()
    Enuc = auto()
    S = auto()
    T = auto()
    V = auto()
    Hcore = auto()
    MM = auto()
    Diff = auto()
    LM = auto()
    AM = auto()
    Eri = auto()

    def __repr__(self):
        return self.name

op_names = list(OperatorNames.__members__.items())

class BaseOperator(np.ndarray):
    """Base operator class.

    Attributes
    ----------
    name : str
        Name of the operator.

    Methods
    -------
    __new__(cls, integral, name)
        Generate new operator.
    
    __init__(self, integral, name)
        Initialize the operator.
    
    assign_name(self,name)
        Assign name to the operator.
    """
    def __new__(cls, integral, name):
        """Generate new operator.

        Parameters
        ----------
        integral : ndarray
            Integral value for the operator.
        
        name : str
            Name of the operator.
        """
        obj = np.asarray(integral).view(cls)
        return obj

    def __init__(self, integral, name):
        """Initialize the operator.

        Parameters
        ----------
        integral : ndarray
            Integral value for the operator.
        
        name : str
            Name of the operator.
        """
        self.assign_name(name)
    
    def assign_name(self,name):
        """Assign name to the operator.

        Parameters
        ----------
        name : OperatorNames
            Name of the operator.

        Raises
        ------
        TypeError
            If name of operator is not OperatorNames.
        """
        if not isinstance(name, OperatorNames):
            raise TypeError("Name of oeprator must be OperatorNames")
        self.name = name

    @classmethod
    def zeros(cls, shape, dtype=float):
        """generate instance with zeros elements ndarray.

        Parameters
        ----------
        shape : tuple/list
            Shape of the ndarray.

        Raises
        ------
        TypeError
            If shape is not tuple or list.
        """
        if not isinstance(shape, tuple or list):
            raise TypeError("Shape of the ndarray must be tuple or list")
        if not len(set(a))==1:
            raise ValueError("Shape of the ndarray must be square/hypersquare")
        obj = np.asarray(np.zeros(shape, dtype=dtype)).view(cls)
        return obj

    @classmethod
    def ones(cls, shape, dtype=float):
        """generate instance with ones elements ndarray.

        Parameters
        ----------
        shape : tuple/list
            Shape of the ndarray.

        Raises
        ------
        TypeError
            If shape is not tuple or list.
        """
        if not isinstance(shape, tuple or list):
            raise TypeError("Shape of the ndarray must be tuple or list")
        if not len(set(a))==1:
            raise ValueError("Shape of the ndarray must be square/hypersquare")
        obj = np.asarray(np.ones(shape, dtype=dtype)).view(cls)
        return obj

    @classmethod
    def random(cls, shape, dtype=float):
        """generate instance with random elements ndarray.

        Parameters
        ----------
        shape : tuple/list
            Shape of the ndarray.

        Raises
        ------
        TypeError
            If shape is not tuple or list.
        """
        if not isinstance(shape, tuple or list):
            raise TypeError("Shape of the ndarray must be tuple or list")
        if not len(set(a))==1:
            raise ValueError("Shape of the ndarray must be square/hypersquare")
        if dtype == float:
            obj = np.random.random(shape)
        elif dtype == complex:
            obj = np.random.random(shape) + np.random.random(shape) * 1j
        else:
            return NotImplementedError('dtype %r not supported!' % dtype)
        obj = np.asarray(obj).view(cls)
        return obj

