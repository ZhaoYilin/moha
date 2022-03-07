from moha.posthf.ci.slater import SlaterDeterminant, SlaterCondon
from moha.posthf.ci.ci_basis_set import CIBasisSet

import numpy as np
import copy
import itertools

__all__ = ['CIHamiltonian']

class CIHamiltonian(np.ndarray):
    """Configuration interaction Hamiltonian.

    Methods
    -------
    __new__(self, integral)
        Generate new CI Hamiltonian.
    
    __init__(self, integral)
        Initialize the CI Hamiltonian.

    Class Methods
    -------------
    build(cls, h1e, g2e, basis_set)
        Build the configuration interaction Hamiltonian.
    
    """
    def __new__(cls, integral):
        """Generate new CI Hamiltonian.

        Parameters
        ----------
        integral : ndarray
            Integral value for the CI Hamiltonian.
        
        Returns
        -------
        ham : CIHamiltonian
            Instance of CIHamiltonian.
        """
        obj = np.asarray(integral).view(cls)
        return obj

    def __init__(self, integral):
        """Initialize the CI Hamiltonian.

        Parameters
        ----------
        integral : ndarray
            Integral value for the CI Hamiltonian.
        """
        pass
    
    @classmethod
    def build(cls, h1e, g2e, basis_set):
        """Build the configuration interaction Hamiltonian.
    
        Parameters
        ----------
        h1e : np.ndarray
            One electron integral.

        g2e : np.ndarray
            Two electron integral.

        basis_set
            Configuration interaction basis set.
        
        Raises
        ------
        TypeError
            If basis_set is not a CIBasisSet instance.

        
        Returns
        -------
        ham : CIHamiltonian
            Instance of CIHamiltonian.
        """
        #Check the input parameters
        if not isinstance(h1e,np.ndarray):
            raise TypeError("Parameter h1e must be a np.ndarray instance.")
        if not isinstance(g2e,np.ndarray):
            raise TypeError("Parameter g2e must be a np.ndarray instance.")
        if not isinstance(basis_set,CIBasisSet):
            raise TypeError("Parameter ci_basis_set must be a CIBasisSet instance.")
        
        #Calculate the ci matrix
        sc = SlaterCondon(h1e,g2e)
        matrix = np.zeros((basis_set.size,basis_set.size))
        for i in range(basis_set.size):
            for j in range(i+1):
                matrix[i,j] = sc.calculate_matrix_element(basis_set[i],basis_set[j])
                matrix[j,i] = matrix[i,j]
        
        #Inital the CIHamiltonian instance
        ham = cls(matrix)
        return ham

