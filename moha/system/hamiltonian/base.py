from moha.system.hamiltonian.abc import ABCHamiltonian
from moha.system.operator.base import BaseOperator
import numpy as np

__all__ = ['BaseHamiltonian']

class BaseHamiltonian(ABCHamiltonian):
    """Base Hamiltonian for a Schrodinger equation.

    Attributes
    ----------
    nspatial : int
        Number of spatial orbitals.

    operators : dict
        Dictionary of operators for the Hamiltonian

    Property
    -------
    npsin : int
        Number of spin orbitals

    Methods
    -------
    __init__(self, nspatial,operators={})
        Initialize the Hamiltonian.

    assign_nspatial(self, nspatial)
        Assigns number of spatial orbitals to the Hamiltonian.

    assign_operators(self, operator)
        Assigns operator integral to the Hamiltonian.

    basis_transformation(self,U)
        Transform the basis set for all the operator integral
    """

    def __init__(self, nspatial,operators={}):
        """Initialize the Hamiltonian.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbitals.

        operators : dict
            Dictionary of operators for the Hamiltonian
        """
        self.assign_nspatial(nspatial)
        self.operators = operators

    @property
    def nspin(self):
        """Return the number of spin orbitals.

        Returns
        -------
        nspin : int
            Number of spin orbitals.

        """
        return self.nspatial*2

    def assign_nspatial(self, nspatial):
        """Assign the number of spatial orbitals.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbital.

        Raises
        ------
        TypeError
            If number of electrons is not an integer.
        ValueError
            If number of electrons is negative number.
        """        
        if not isinstance(nspatial, int):
            raise TypeError("Number of spatial orbitals must be an integer")
        if nspatial < 0:
            raise ValueError("Number of electron must be a none negative integer")
        self.nspatial = nspatial

    def assign_operator(self, operator):
        """Assign operator.

        Parameters
        ----------
        operator : Operator instance
            Operator for the Hamiltonain

        Raises
        ------
        TypeError
            If number of electrons is not an integer.
        """
        if not isinstance(operator, BaseOperator):
            raise TypeError("operator must be a BaseOperator instance")
        self.operators[operator.name] = operator

    def basis_transformation(self,matrix):
        """Transform the basis set for all the operator integral

        Parameters
        ----------
        matrix : np.ndarray(K, K)
            Transformation matrix
        """
        for name,operator in self.operators.items():
            if isinstance(operator,(OneElectronOperator,TwoElectronOperator)):
                operator.basis_transformation(matrix)
