import numpy as np
class PopulationAnalysis(object):
    """
    """
    def __init__(self,charges=[]):
        self.charges = charges

    @classmethod
    def mulliken(cls,molecule,basis_set,wavefunction,hamiltonian):
        """
        Mulliken population analysis
        D: density matrix
        S: overlap matrix
        """
        population = cls()
        S = hamiltonian.operators['overlap'].value
        D = wavefunction.density
        if type(D) is np.ndarray:
            DS = np.dot(D,S)
        elif type(D) is dict:
            DS = np.dot(D['alpha'],S) + np.dot(D['beta'],S)
        for i,atom in enumerate(molecule.atoms):
            charge = 0.
            for j,orbital in enumerate(basis_set.basis):
                if orbital.atom_index==i:
                    charge += DS[j,j]
            charge = atom.number - charge
            print('{0:2s} {1:3d} {2:4s} {3:5f}'.format('Charge on atom', i, ':', charge))
            population.charges.append(charge)
        return population

    @classmethod
    def lowdin(cls,molecule,basis_set,wavefunction,hamiltonian):
        """
        Lowdin population analysis
        """
        population = cls()
        D = wavefunction.density
        S = hamiltonian.operators['overlap'].value
        val, vec = np.linalg.eig(S)
        val_minus_half = (np.diag(val**(0.5)))
        X = np.dot(vec,np.dot(val_minus_half,np.transpose(vec)))
        if type(D) is np.ndarray:
            XDX = np.dot(X,np.dot(D,X))
        elif type(D) is dict:
            XDX = np.dot(X,np.dot(D['alpha'],X))+np.dot(X,np.dot(D['beta'],X))
        for i,atom in enumerate(molecule.atoms):
            charge = 0.
            for j,orbital in enumerate(basis_set.basis):
                if orbital.atom_index==i:
                    charge += XDX[j,j]
            charge = atom.number - charge
            print('{0:2s} {1:3d} {2:4s} {3:5f}'.format('Charge on atom', i, ':', charge))
            population.charges.append(charge)
        return population
