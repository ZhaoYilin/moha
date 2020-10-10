from abc import ABC,abstractmethod,abstractproperty
import numpy as np
from moha.system.operator import *
from moha.log import timer

__all__ = ['ABCHamiltonian']

class ABCHamiltonian(ABC):
    """Abstract Hamiltonian for a Schrodinger equation.

    Attributes
    ----------
    nspatial : int
        Number of spatial orbitals.

    operators : dict
        Dictionary of operators for the Hamiltonian

    Abstract Property
    -----------------
    npsin : int
        Number of spin orbitals

    Methods
    -------
    __init__(self, nspatial,operators={})
        Initialize the Hamiltonian.

    Abstract Methods
    ----------------
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

    @abstractproperty
    def nspin(self):
        """Return the number of spin orbitals.

        Returns
        -------
        nspin : int
            Number of spin orbitals.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError

    @abstractmethod
    def assign_nspatial(self, nspatial):
        """Assign the number of spatial orbitals.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbital.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError

    @abstractmethod
    def assign_operator(self, operator):
        """Assign operator.

        Parameters
        ----------
        operator : Operator instance
            Operator for the Hamiltonain

        Raises
        ------
        NotImplementedError
            If called.

        """
        raise NotImplementedError

    @abstractmethod
    def basis_transformation(self,U):
        """Transform the basis set for all the operator integral

        Parameters
        ----------
        U : np.ndarray(K, K)
            Transformation matrix

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError

