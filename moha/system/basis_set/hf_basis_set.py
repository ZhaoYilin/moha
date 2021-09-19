from moha.system.basis_set.base import BaseBasisSet
from moha.system.basis.gaussian_orbital import GaussianOrbital

__all__ = ['HFBasisSet']

class HFBasisSet(BaseBasisSet):
    """Hatree Fock basis set.
    
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
        super().__init__(size,bases)

    def append_basis(self,basis):
        """Add one basis to the basis set.

        Parameters
        ----------
        basis
            Basis instance.
        
        Raises
        ------
        TypeError
            If orb is not a GaussianOrbital instance.
        """    
        if not isinstance(basis, GaussianOrbital):
            raise TypeError("Basis must be a GaussianOrbital instance")
        self.size += 1
        self.bases.append(basis)
