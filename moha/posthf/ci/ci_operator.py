from moha.posthf.ci.slater import SlaterDeterminant, SlaterCondon
from moha.posthf.ci.ci_basis_set import CIBasisSet
from moha.system.operator.base import OperatorNames, BaseOperator

import numpy as np
import copy

__all__ = ['CIOperator']

class CIOperator(BaseOperator):
    """Configuration interaction Operator.

    Methods
    -------
    __new__(self, integral)
        Generate new CI Operator.
    
    __init__(self, integral)
        Initialize the CI Operator.

    Class Methods
    -------------
    build(cls, h1e, g2e, basis_set)
        Build the CI Operator.
    
    """
    def __init__(cls, integral, name):
        """Generate new CI Operator.

        Parameters
        ----------
        integral : ndarray
            Integral value for the CI Operator.
        
        Returns
        -------
        ham : CIOperator
            Instance of CIOperator.
        """
        super().__init__(integral, name)

    @classmethod
    def build(cls, h1e, g2e, basis_set):
        """Build the configuration interaction Operator.
    
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
        ham : CIOperator
            Instance of CIOperator.
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
        integral = np.zeros((basis_set.size,basis_set.size))
        for i in range(basis_set.size):
            for j in range(i+1):
                integral[i,j] = sc.calculate_matrix_element(basis_set[i],basis_set[j])
                integral[j,i] = integral[i,j]
        
        #Inital the CIOperator instance
        operator = cls(integral, OperatorNames.CI)
        return operator

