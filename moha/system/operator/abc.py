from abc import ABC,abstractmethod,abstractproperty

class ABCOperator(ABC):
    """Abstract operator class.

    Attributes
    ----------
    name : str
        Name of the operator.

    nspatial : int
        Number of the spatial orbitals.

    nelectron : int    
        Number of electrons.

    integral
        Integral of the operator.

    Abstract Property
    -----------------
    nspin : int
        Number of spin orbitals.
    
    dtype 
        Data type of the integrals.

    Methods
    -------
    __init__(self,name,nspatial,nelectron=None,integral=None)
        Initialize the Operator.        
    
    Abstract Method
    ---------------
    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.
    
    assign_nelectron(self,nelectron)
        Assign number of elctrons involved in the integral.

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
        
    @abstractproperty
    def nspin(self):
        """Number of spin orbitals.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError

    @abstractproperty
    def dtype(self):
        """Data type of the integrals.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError

    @abstractmethod
    def assign_name(self,name):
        """Assign name to the operator.

        Parameters
        ----------
        name : str
            Name of the operator.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError 

    @abstractmethod
    def assign_nspatial(self,nspatial):
        """Assign number of spatial orbitals for the operator.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbitals for the operator.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError

    @abstractmethod
    def assign_nelectron(self,nelectron):
        """Assign number of electrons involved in the integral.

        Parameters
        ----------
        nelectron : int
            Number of electrons involved in the integral.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError

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

