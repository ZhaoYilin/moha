import numpy as np
from moha.system.operator import *
from moha.system.hamiltonian.base import BaseHamiltonian
from moha.io.log import timer

__all__ = ['ChemicalHamiltonian']

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
    
    build(cls,molecule,basis_set)
        Build the Chemical Hamiltonian

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
        super().__init__(nspatial,operators)

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
        nuclear_repulsion = NuclearRepulsionOperator.build(molecule,basis_set)
        overlap = OverlapOperator.build(basis_set)
        kinetic = KineticOperator.build(basis_set)
        nuclear_attraction = NuclearAttractionOperator.build(molecule,basis_set)
        electron_repulsion = ElectronRepulsionOperator.build(basis_set)
        
        h1e_integral = kinetic.integral + nuclear_attraction.integral
        g2e_integral = electron_repulsion.integral
        h1e = OneElectronOperator('h1e',nspatial,nelectron=1,integral=h1e_integral)
        g2e = TwoElectronOperator('g2e',nspatial,nelectron=2,integral=g2e_integral)

        #assign operators to chemical hamiltonian
        ham.assign_operator(nuclear_repulsion)
        ham.assign_operator(overlap)
        ham.assign_operator(kinetic)
        ham.assign_operator(nuclear_attraction)
        ham.assign_operator(electron_repulsion)
        
        ham.assign_operator(h1e)
        ham.assign_operator(g2e)

        return ham
    
    @classmethod
    @timer.with_section('Integral')
    def from_fcidump(cls,fcidump):
        """Build the chemical hamiltonian

        Parmeters
        ---------
        molecule : Molecule instance
            Instance of molecule class

        basis_set : Basis Set instance
            Instance of one electron basis set
        """
        nspatial = fcidump.n_sites
        ham = cls(nspatial)
        
        #build operators
        nuclear_repulsion = NuclearRepulsionOperator('nuclear_repulsion',nspatial,nelectron=0,integral=fcidump.const_e)
        h1e = OneElectronOperator('h1e',nspatial,nelectron=1,integral=fcidump.h1e)
        g2e = TwoElectronOperator('g2e',nspatial,nelectron=2,integral=fcidump.g2e)
        
        #assign operators to chemical hamiltonian
        ham.assign_operator(nuclear_repulsion)
        ham.assign_operator(h1e)
        ham.assign_operator(g2e)
        
        return ham

