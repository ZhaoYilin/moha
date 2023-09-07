from moha.io.log import log,timer

import numpy as np

__all__ = ['Population']

class Population(object):
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

    @timer.with_section('Population')
    def mulliken(self):
        """Kernel of the solver.
        
        Returns
        -------
        results : dict
            Mulliken population analysis results.
        """
        log.hline(char='=')
        log.blank()
        log('Population Analysis Section'.format())

        S = self.ham['S']
        D = self.wfn.density_matrix
        DS = np.dot(D,S)

        population = []

        for i,atom in enumerate(self.mol):
            charge = 0.
            for j,basis in enumerate(self.basis_set):
                if basis.origin==atom.coordinate:
                    charge += DS[j,j]
            charge = atom.number - 2*charge
            population.append(charge)

        log.hline()
        log('Mulliken Population'.format())
        log.hline()
        for i,charge in enumerate(population):
            log('{0:<15s} {1:<4d} {2:<20.6f}'.format('Charge on atom', i, charge))
        log('{0:<20s} {1:<20.6f}'.format('Total Charge', sum(population)))
        log.blank()

        results = {
        "success": True,
        "population": population
        }
        return results

    @timer.with_section('Population')
    def lowdin(self):
        """Kernel of the solver.
        
        Returns
        -------
        results : dict
            Mulliken population analysis results.
        """
        log.hline(char='=')
        log.blank()
        log('Population Analysis Section'.format())

        S = self.ham['S']
        D = self.wfn.density_matrix
        val, vec = np.linalg.eig(S)
        val_minus_half = (np.diag(val**(0.5)))
        X = np.dot(vec,np.dot(val_minus_half,np.transpose(vec)))
        XDX = np.dot(X,np.dot(D,X))

        population = []

        for i,atom in enumerate(self.mol):
            charge = 0.
            for j,basis in enumerate(self.basis_set):
                if basis.origin==atom.coordinate:
                    charge += XDX[j,j]
            charge = atom.number - 2*charge
            population.append(charge)

        log.hline()
        log('Lowdin Population'.format())
        log.hline()
        for i,charge in enumerate(population):
            log('{0:<15s} {1:<4d} {2:<20.6f}'.format('Charge on atom', i, charge))
        log('{0:<20s} {1:<20.6f}'.format('Total Charge', sum(population)))
        log.blank()

        results = {
        "success": True,
        "population": population
        }

        return results

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
