from abc import ABC,abstractmethod,abstractproperty
import numpy as np
from moha.system.operator import *
from moha.io.log import timer

__all__ = ['ChemicalHamiltonian']

class BaseHamiltonian(ABC):
    """Base Hamiltonian for a Schrodinger equation.

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

class ChemicalHamiltonian(BaseHamiltonian):
    """Chemical Hamiltonian for a Schrodinger equation.

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

    @classmethod
    @timer.with_section('Integral')
    def build(cls,molecule,basis_set):
        """Build the chemical hamiltonian

        Parmeters
        ---------
        molecule : Molecule instance
            Instance of molecule class

        basis_set : Basis Set instance
            Instance of one electron basis set
        """
        nspatial = basis_set.size
        ham = cls(nspatial)
        #build operators
        nuclear_repulsion = ZeroElectronOperator.build('nuclear_repulsion',molecule)
        overlap = OneElectronOperator.build('overlap',molecule,basis_set)
        kinetic = OneElectronOperator.build('kinetic',molecule,basis_set)
        nuclear_attraction = OneElectronOperator.build('nuclear_attraction',molecule,basis_set)
        electron_repulsion = TwoElectronOperator.build('electron_repulsion',molecule,basis_set)
        #assign operators to chemical hamiltonian
        ham.assign_operator(nuclear_repulsion)
        ham.assign_operator(overlap)
        ham.assign_operator(kinetic)
        ham.assign_operator(nuclear_attraction)
        ham.assign_operator(electron_repulsion)

        return ham

    @property
    def nspin(self):
        """Return the number of spin orbitals.

        Returns
        -------
        nspin : int
            Number of spin orbitals.

        """
        self.nspin = nspatial*2

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
            If number of electrons is not a positive number.
        """        
        if not isinstance(nspatial, int):
            raise TypeError("Number of spatial orbitals must be an integer")
        if nspatial <= 0:
            raise ValueError("Number of spatial orbitals must be a positive integer")
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

    def basis_transformation(self,U):
        """Transform the basis set for all the operator integral

        Parameters
        ----------
        U : np.ndarray(K, K)
            Transformation matrix
        """
        for name,operator in self.operators.items():
            if isinstance(operator,(OneElectronOperator,TwoElectronOperator)):
                operator.basis_transformation(U)
