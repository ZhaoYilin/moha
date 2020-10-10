from abc import abstractmethod,abstractproperty
from moha.system.operator.abc import ABCOperator
import numpy as np

class BaseOperator(ABCOperator):
    """Base operator class.

    Attributes
    ----------
    name : str
        Name of the operator.

    nspatial : int
        Number of the spatial orbitals.
    
    nelectron : int
        Number of electrons.

    integral : int,float,np.ndarray,2-dictionary of np.ndarray
        Integral of the operator.

    Property
    --------
    npsin : int
        Number of spin orbitals.
    
    Abstract Property
    -----------------
    dtype
        Data type of the integrals.

    Methods
    -------
    __init__(self,name,nspatial,nelectron=None,integral=None)
        Initialize the Operator.        
    
    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.
    
    assign_nelectron(self,nelectron)
        Assign number of elctrons involved in the integral.

    Abstract Methods
    ----------------
    assign_integral(self,integral)
        Assign integral for the operator.
    """
    def __init__(self,name,nspatial,nelectron=None,integral=None):
        """Initialize the operator

        Parameters
        ----------
        name : str
            Name of the operator.

        nspatial : int
            Number of the spatial orbitals for the operator

        nelectron : int
            Number of elctrons involved in the integral.

        integral
            Integral value for the operator.

        """
        self.assign_name(name)
        self.assign_nspatial(nspatial)
        self.assign_nelectron(nelectron)
        self.assign_integral(integral)

    @property        
    def nspin(self):
        """Return the number of spin orbitals.

        Returns
        -------
        nspin : int
            Number of spin orbitals.
        """
        return self.nspatial*2

    @abstractproperty
    def dtype(self):
        """Data type of the integrals.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError

    def assign_name(self,name):
        """Assign name to the operator.

        Parameters
        ----------
        name : str
            Name of the operator.

        Raises
        ------
        TypeError
            If name of operator is not str.
        """
        if not isinstance(name, str):
            raise TypeError("Name of oeprator must be str")
        self.name = name

    def assign_nspatial(self,nspatial):
        """Assign number of spatial orbitals for the operator.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbitals for the operator.

        Raises
        ------
        TypeError
            If number of spatial orbitals is not an integer.

        ValueError
            If number of spatial orbitals is not a positive number.
        """
        if not isinstance(nspatial, int):
            raise TypeError("Number of spatial orbitals must be an integer")
        if nspatial <= 0:
            raise ValueError("Number of spatial orbitals must be a positive number")
        self.nspatial = nspatial
    
    def assign_nelectron(self,nelectron):
        """Assign number of elctrons involved in the integral.

        Parameters
        ----------
        nelec : int
            Number of elctrons involved in the integral.

        Raises
        ------
        TypeError
            If number of electrons is not an integer.

        ValueError
            If number of electrons is not a positive number.
        """
        if not isinstance(nelectron, int):
            raise TypeError("Number of electrons must be an integer")
        if nelectron < 0:
            raise ValueError("Number of electrons must be a none negative number")
        self.nelectron = nelectron

    @abstractmethod
    def assign_integral(self,integral):
        """Assign integral for the operator.

        Parameters
        ----------
        integral
            Integral for the operator.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError
