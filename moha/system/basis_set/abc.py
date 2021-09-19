from abc import ABC,abstractmethod,abstractproperty

__all__ = ['ABCBasisSet']

class ABCBasisSet(ABC):
    """Abstract basis set.
    
    Attributes
    ----------
    size : int
        Size of the basis set.

    bases : list
        A list of basis.

    Methods
    -------
    __init__(self,size=0,bases=[])
        Initialize the basis set.

    Abstract Methods
    ----------------
    assign_size(self,size)
        Assign size to the basis set.
    
    assign_bases(self,bases)
        Assign bases to the basis set.
    
    append_basis(self,basis)
        Add one basis tp the basis set.
    """
    def __init__(self,size=0,bases=[]):
        """Initialize the instance.

        Parameters
        ----------
        size : int
            Size of the basis set.

        bases : list
            A list of basis.
        """ 
        self.assign_size(size)
        self.assign_bases(bases)
    
    @abstractmethod
    def assign_size(self,size):
        """Assign size to the basis set.

        Parameters
        ----------
        size : int
            Size of the basis set.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError        
    
    @abstractmethod
    def assign_bases(self,bases):
        """Assign bases to the basis set.

        Parameters
        ----------
        bases : list
            A list of basis.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError        

    @abstractmethod
    def append_basis(self,basis):
        """Add one basis to the basis set.

        Parameters
        ----------
        basis
            Basis instance.
        
        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError        
