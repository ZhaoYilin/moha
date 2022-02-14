from moha.symmetry.symmetry import SZ
from moha.system.operator import *
from moha.system.hamiltonian.base import BaseHamiltonian
from moha.io.iofcidump import FCIDUMP
from moha.io.log import timer

import numpy as np

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
        Enuc = NuclearRepulsionOperator.build(molecule)
        S = OverlapOperator.build(basis_set)
        T = KineticOperator.build(basis_set)
        V = NuclearAttractionOperator.build(molecule,basis_set)
        Hcore = OneElectronOperator(T+V,OperatorNames.Hcore)
        Eri = ElectronRepulsionOperator.build(basis_set)
        
        #assign operators to chemical hamiltonian
        ham.assign_operator(Enuc)
        ham.assign_operator(S)
        ham.assign_operator(T)
        ham.assign_operator(V)
        ham.assign_operator(Hcore)
        ham.assign_operator(Eri)

        return ham
    
    @classmethod
    @timer.with_section('Integral')
    def from_numpy(cls,Enuc,S,h1e,g2e):
        """Build the chemical hamiltonian

        Parmeters
        ---------
        molecule : Molecule instance
            Instance of molecule class

        basis_set : Basis Set instance
            Instance of one electron basis set
        """
        nspatial = h1e.shape[0]
        ham = cls(nspatial)
        
        #build operators
        Enuc = NuclearRepulsionOperator(Enuc,OperatorNames.Enuc)
        S = OneElectronOperator(S,OperatorNames.S)
        Hcore = OneElectronOperator(h1e,OperatorNames.Hcore)
        Eri= ElectronRepulsionOperator(g2e,OperatorNames.Eri)

        #assign operators to chemical hamiltonian
        ham.assign_operator(Enuc)
        ham.assign_operator(S)
        ham.assign_operator(Hcore)
        ham.assign_operator(Eri)
        
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
        self.fcidump = fcidump
        self.orb_sym = fcidump.orb_sym
        self.n_syms = max(self.orb_sym) + 1
        self.nspatial = fcidump.n_sites
        self.vacuum = SZ(0, 0, 0)
        self.target = SZ(fcidump.n_elec, fcidump.twos, fcidump.ipg)
        
        ham = cls(self.nspatial)
        
        #build operators
        Enuc = NuclearRepulsionOperator(fcidump.const_e,OperatorNames.Enuc)
        Hcore = OneElectronOperator(fcidump.h1e,OperatorNames.Hcore)
        Eri= ElectronRepulsionOperator(fcidump.g2e,OperatorNames.Eri)

        #assign operators to chemical hamiltonian
        ham.assign_operator(Enuc)
        ham.assign_operator(Hcore)
        ham.assign_operator(Eri)
        
        return ham

