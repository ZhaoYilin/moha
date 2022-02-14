from moha.io.log import log,timer
from moha.system.operator.base import OperatorNames

import numpy as np

__all__ = ['PopulationAnalysisMulliken','PopulationAnalysisLowdin']

class PopulationAnalysis(object):
    """Population Analysis solver.

    Attributes
    ----------
    mol
        Chemical molecule.
    
    basis_set
        Basis set of the molecule
    
    ham
        Chemical Hamiltonian.

    wfn
        Hartree Fock wavefunction.

    Methods
    -------
    __init__(self,mol,basis_set,ham,wfn)
        Initialize the solver.
    
    assign_molecule(self,mol)
        Assign chemical molecule to the solver.
    
    assign_basis_set(self,basis_set)
        Assign basis set to the solver.

    assign_hamiltonian(self,ham)
        Assign chemical Hamiltonian to the solver.

    assign_wavefunction(self,wfn)
        Assign Hartree Fock wavefunction to the solver.
    """
    def __init__(self,mol,basis_set,ham,wfn):
        """Initialize the solver.

        Attributes
        ----------
        mol
            Chemical molecule.
        
        basis_set
            Basis set of the molecule
            
        ham
            Chemical Hamiltonian.

        wfn
            Hartree Fock wavefunction.
        """
        self.assign_molecule(mol)
        self.assign_basis_set(basis_set)
        self.assign_hamiltonian(ham)
        self.assign_wavefunction(wfn)
    
    def assign_molecule(self,mol):
        """Assign chemical molecule to the solver.

        Attributes
        ----------
        mol
            Chemical molecule.
        """
        self.mol = mol
    
    def assign_basis_set(self,basis_set):
        """Assign basis set to the solver.

        Attributes
        ----------
        basis_set
            Basis set of the molecule.
        """
        self.basis_set = basis_set

    def assign_hamiltonian(self,ham):
        """Assign chemical Hamiltonian to the solver.

        Attributes
        ----------
        ham
            Chemical Hamiltonian.
        """
        self.ham = ham

    def assign_wavefunction(self,wfn):
        """Assign Hartree Fock wavefunction to the solver.

        Attributes
        ----------
        wfn
            Hartree Fock wavefunction.
        """
        self.wfn = wfn

class PopulationAnalysisMulliken(PopulationAnalysis):
    """Mulliken population analysis.
    
    Attributes
    ----------
    mol
        Chemical molecule.
    
    basis_set
        Basis set of the molecule
    
    ham
        Chemical Hamiltonian.

    wfn
        Hartree Fock wavefunction.

    Methods
    -------
    __init__(self,mol,basis_set,ham,wfn)
        Initialize the solver.
    
    assign_molecule(self,mol)
        Assign chemical molecule to the solver.
    
    assign_basis_set(self,basis_set)
        Assign basis set to the solver.

    assign_hamiltonian(self,ham)
        Assign chemical Hamiltonian to the solver.

    assign_wavefunction(self,wfn)
        Assign Hartree Fock wavefunction to the solver.
    """
    def __init__(self,mol,basis_set,ham,wfn):
        super().__init__(mol,basis_set,ham,wfn)

    @timer.with_section('Population')
    def kernel(self):
        """Kernel of the solver.
        
        Returns
        -------
        results : dict
            Mulliken population analysis results.
        """
        log.hline()
        log('Population Analysis Section'.format())
        log.hline()

        S = self.ham.operators[OperatorNames.S]
        D = self.wfn.density_matrix
        DS = np.dot(D,S)
        
        population = []

        for i,atom in enumerate(self.mol.atoms):
            charge = 0.
            for j,basis in enumerate(self.basis_set.bases):
                if basis.atom_index==i:
                    charge += DS[j,j]
            charge = atom.number - charge
            population.append(charge)

        log.hline()
        log('Mulliken Population'.format())
        for i,charge in enumerate(population):
            log('{0:2s} {1:3d} {2:4s} {3:5f}'.format('Charge on atom', i, ':', charge))
        log.hline()
        
        results = {
        "success": True,
        "population": population
        }
        return results

class PopulationAnalysisLowdin(PopulationAnalysis):
    """Lowdin population analysis.
    
    Attributes
    ----------
    mol
        Chemical molecule.
    
    basis_set
        Basis set of the molecule
    
    ham
        Chemical Hamiltonian.

    wfn
        Hartree Fock wavefunction.

    Methods
    -------
    __init__(self,mol,basis_set,ham,wfn)
        Initialize the solver.
    
    assign_molecule(self,mol)
        Assign chemical molecule to the solver.
    
    assign_basis_set(self,basis_set)
        Assign basis set to the solver.

    assign_hamiltonian(self,ham)
        Assign chemical Hamiltonian to the solver.

    assign_wavefunction(self,wfn)
        Assign Hartree Fock wavefunction to the solver.
    """
    def __init__(self,mol,basis_set,ham,wfn):
        super().__init__(mol,basis_set,ham,wfn)

    @timer.with_section('Population')
    def kernel(self):
        """Kernel of the solver.
        
        Returns
        -------
        results : dict
            Mulliken population analysis results.
        """
        log.hline()
        log('Population Analysis Section'.format())
        log.hline()

        S = self.ham.operators[OperatorNames.S]
        D = self.wfn.density_matrix
        val, vec = np.linalg.eig(S)
        val_minus_half = (np.diag(val**(0.5)))
        X = np.dot(vec,np.dot(val_minus_half,np.transpose(vec)))
        XDX = np.dot(X,np.dot(D,X))
        
        population = []
        
        for i,atom in enumerate(self.mol.atoms):
            charge = 0.
            for j,basis in enumerate(self.basis_set.bases):
                if basis.atom_index==i:
                    charge += XDX[j,j]
            charge = atom.number - charge
            population.append(charge)

        log.hline()
        log('Lowdin Population'.format())
        for i,charge in enumerate(population):
            log('{0:2s} {1:3d} {2:4s} {3:5f}'.format('Charge on atom', i, ':', charge))
        log.hline()
        
        results = {
        "success": True,
        "population": population
        }
        return results

