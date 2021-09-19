from moha.system.basis_set.abc import ABCBasisSet

__all__ = ['BaseBasisSet']

class BaseBasisSet(ABCBasisSet):
    """Base basis set.
    
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
    
    def assign_size(self,size):
        """Assign size to the basis set.

        Parameters
        ----------
        size : int
            Size of the basis set.

        Raises
        ------
        TypeError
            If size of the basis set is not an integer.

        ValueError
            If size of the basis set is not a positive number.        
        """
        if not isinstance(size, int):
            raise TypeError("Size of the basis set must be an integer")
        if size < 0:
            raise ValueError("Size of the basis set must be a none negative number")
        self.size = size
    
    def assign_bases(self,bases):
        """Assign bases to the basis set.

        Parameters
        ----------
        bases : list
            A list of basis.
        
        Raises
        ------
        TypeError
            If bases is not a list.
        """
        if not isinstance(bases, list):
            raise TypeError("Bases must be a list")
        self.bases = bases

    def append_basis(self,basis):
        """Add one basis to the basis set.

        Parameters
        ----------
        basis
            Basis instance.
        """    
        self.size += 1
        self.bases.append(basis)
