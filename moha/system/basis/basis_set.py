from moha.system.basis.gaussian_orbital import GaussianOrbital

__all__ = ['BasisSet']

class BasisSet(list):
    """Hatree Fock basis set.
    
    Attributes
    ----------

    Methods
    -------
    __new__(cls)
        Generate new basis set.
    
    __init__(self)
        Initialize the basis set.
    
    append(self,basis)
        Add basis to the basis set.
    """
    def __new__(cls):
        """Generate new basis set.
        """
        obj = list.__new__(cls)
        
        return obj

    def __init__(self):
        """Initialize the basis set.
        """
        pass
    
    def append(self,basis):
        """Add basis to the basis set.

        Parameters
        ----------
        basis : GaussianOrbital
            Gaussian type orbital.
        
        Raises
        ------
        TypeError
            If orb is not a GaussianOrbital instance.
        """
        if not isinstance(basis, GaussianOrbital):
            raise TypeError("Basis must be a GaussianOrbital instance")
        super().append(basis)

