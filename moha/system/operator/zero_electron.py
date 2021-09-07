from moha.system.operator.base import BaseOperator
from moha.system.molecule import Molecule
from moha.system.basis_set import GeneralBasisSet
import numpy as np

class NuclearRepulsionOperator(BaseOperator):
    """Nuclear repulsion operator class.

    Attributes
    ----------
    name : str
        Name of the operator.

    nspatial : int
        Number of the spatial orbitals.

    nelectron : int
        Number of electrons.

    integral : int,float
        Integral of the operator.

    Property
    --------
    npsin : int
        Number of spin orbitals.

    dtype : int,float
        Data type of the integrals.

    Methods
    -------
    __init__(self,name,nspatial,nelectron,integral)
        Initialize the Operator.

    build(cls,molecule,basis_set)
        Build the operator instance.

    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.

    assign_nelectron(self,nelectron)
        Assign number of elctrons involved in the integral.

    assign_integral(self,integral)
        Assign integral for the operator.
    """

    def __init__(self,name,nspatial,nelectron=0,integral=None):
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
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,molecule,basis_set):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        moleucle : Molecule 
            Molecule instance.

        basis_set : GeneralBasisSet
            General basis set instance.

        Raises
        ------
        TypeError
            If molecule parameter is not a Molecule instance.
            If basis_set parameter is not a GeneralBasisSet instance.
        """
        if not isinstance(molecule,Molecule):
            raise TypeError("molecule parameter must be a Molecule instance.")        
        if not isinstance(basis_set,GeneralBasisSet):
            raise TypeError("basis_set parameter must be a GeneralBasisSet instance.")        
        nspatial = basis_set.size
        E = 0.
        for i,A in enumerate(molecule.atoms):
            for j,B in enumerate(molecule.atoms):
                if j>i:
                    E += (A.number*B.number)/molecule.bond_length(i,j)
        operator = cls('nuclear_repulsion',nspatial,0,E)

        return operator

    @property
    def dtype(self):
        """Data type of the integrals.

        Returns
        -------
        dtype : int,float
            Data type of the integrals.

        """
        return self.integral.dtype

    def assign_integral(self,integral):
        """Assign integral value for the operator.

        Parameters
        ----------
        integral : int, float
            Integral value for the operator.

        Raises
        ------
        TypeError
            If value is not an integer or float number.
        ValueError
            If value is not none negative.
        """
        if integral is None:
            self.integral = 0.0
        if not isinstance(integral, (int,float)):
            raise TypeError("Integral must be integer or float.")
        if integral < 0:
            raise ValueError("Integral must be none negative.")
        self.integral = integral        
