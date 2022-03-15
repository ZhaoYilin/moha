from moha.system.operator.base import OperatorNames, BaseOperator
from moha.system.molecule import Molecule

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
    def __init__(self, integral, name=OperatorNames.Enuc):
        """Initialize the operator

        Parameters
        ----------
        integral : ndarray
            Integral value for the operator.
        
        name : OperatorNames
            name of the operator.
        """
        super().__init__(integral, name)
    
    @classmethod
    def build(cls,molecule):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        moleucle : Molecule 
            Molecule instance.

        Raises
        ------
        TypeError
            If molecule parameter is not a Molecule instance.
        """
        if not isinstance(molecule,Molecule):
            raise TypeError("molecule parameter must be a Molecule instance.")        
        E = 0.
        for i,A in enumerate(molecule):
            for j,B in enumerate(molecule):
                if j>i:
                    E += (A.number*B.number)/molecule.bond_length(i,j)
        operator = cls(E,OperatorNames.Enuc)

        return operator
