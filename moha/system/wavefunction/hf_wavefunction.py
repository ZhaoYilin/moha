from moha.system.wavefunction.wavefunction import WaveFunction
from moha.system.basis.basis import SlaterDeterminant
import numpy as np 

__all__ = ['HFWaveFunction']

class HFWaveFunction(WaveFunction):
    """Hatree Fock wavefunction class.

    Attributes
    ----------
    nelec : int
        Number of electrons.
    
    nspatial : int
        Number of spatial orbitals.
    
    basis_set : Basis
        Basis set of the wavefunction.
    
    density_matrix : np.ndarray
        Density matrix of the wavefunction.
    
    orbital_energies : np.ndarray
        Orbital energies of the wavefunction.

    Properties
    ----------
    spin : int
        Spin of the wavefunction

    Methods
    -------
    __init__(self, nelec, nspatial, basis_set=None, coefficients=None)
        Initialize the wavefunction.
    
    assign_nelec(self, nelec)
        Assign the number of electrons.
    
    assign_nspatial(self, nspatial)
        Assign the number of spatial orbitals.

    assign_occ(self, occ)
        Assign the occupation number of the wavefunction.    
    
    assign_basis_set(self, basis_set)
        Assign basis set of the wavefunction.
    
    assign_coefficients(self, coefficients)
        Assign coefficients of the wavefunction.
    
    assign_density_matrix(self, density_matrix)
        Assign density matrix of the wavefunction.
    
    assign_orbital_energies(self, orbital_energies)
        Assign orbital energies of the wavefunction.
    """
    @classmethod
    def build(cls,basis_set,symmetry):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : BasisSet
            Hatree Fock basis set instance.

        Raises
        ------
        TypeError
            If basis_set parameter is not a BasisSet instance.

        Return
        ------
        wfn: HFWaveFunction
            Return a hantree fock wavefunction instance.
        """
        norb = len(basis_set)
        ne = symmetry.n
        ms2 = symmetry.ms2
        na = int((ne+ms2)/2)
        nb = int((ne-ms2)/2)

        alpha = [1]*na + [0]*(nspa-na)
        beta = [1]*nb + [0]*(nspatials-nb)

        configuration = []
        for i,j in zip(alpha, beta):
            configuration.append(i)
            configuration.append(j)

        ansatz = SlaterDeterminant(configuration)
        c = np.identity((nspatial,nspatial))

        wfn = cls()
        wfn.assign_ansatz(ansatz)
        wfn.assign_basis_set(basis_set)
        wfn.assign_coefficient(coefficient)
        wfn.assign_symmetry(symmetry)

        return wfn

    def assign_ansatz(self,ansatz):
        """Assign ansaztz to the wave function.

        Parameters
        ----------
        ansatz : SlaterDeterminant
            The SlaterDeterminant instance.
        
        Raises
        ------
        TypeError
            If ansatz is not SlaterDeterminant instance.
        """
        if not isinstance(ansatz,SlaterDeterminant):
            raise TypeError("Parameters ansatz must be SlaterDeterminant instance.")
        self.ansatz = ansatz

    def assign_basis_set(self, basis_set):
        """Assign the basis_set of the wave function.

        Parameters
        ----------
        basis_set
            Basis set of the wave function.

        Raises
        ------
        TypeError
            If basis_set is not BasisSet instance.
        """
        from moha.system.basis.basis_set import BasisSet
        if not isinstance(basis_set,BasisSet):
            raise TypeError("Parameters basis_set must be BasisSet instance.")
        self.basis_set = basis_set

    def assign_coefficient(self,coefficient):
        """Assign coefficient to the wave function.

        Parameters
        ----------
        coefficient: np.ndarray
            Coefficient of the wave function.

        Raises
        ------
        TypeError
            If coefficient is not np.ndarray instance.
        """
        if not isinstance(coefficient,np.ndarray):
            raise TypeError("Parameters coefficient must be np.ndarray instance.")
        self.coefficient = coefficient

    def assign_symmetry(self,symmetry):
        """Assign symmetry to the wave function.

        Parameters
        ----------
        symmetry: Symmetry instance
            Symmetry of the wave function.

        Raises
        ------
        TypeError
            If symmetry is not Symmetry instance.
        """
        from moha.symmetry.symmetry import Symmetry
        if not isinstance(symmetry,Symmetry):
            raise TypeError("Parameters symmetry must be Symmetry instance.")
        self.symmetry = symmetry

    @property
    def nelec(self):
        """The exact state is a function of the space and spin coordinates of N
        electrons.
        """
        return self.symmetry.n

    @property
    def spin(self,ms2):
        """In nonrelativistic theory, the exact stationary states are 
        eigenfunctions of the total and projected spin operators.
        """
        return self.symmetry.ms2

    @property
    def spatial_symmetry(self,ipg):
        """Within the Born-Oppenheimer approximation, the exact stationary 
        states form a basis for an irreducible representation of the 
        molecular point group.
        """
        return self.symmetry.ipg

    def assign_density_matrix(self, density_matrix):
        """ Assign density matrix of the wavefunction.

        Parameters
        ----------
        density_matrix : np.ndarray
            Density matrix of the wavefunction.

        Raises
        ------
        """
        #if density_matrix is None:
        #    density_matrix = np.zeros((self.nspatial,self.nspatial))
        #elif not isinstance(density_matrix,np.ndarray):
        #    raise TypeError("Density matrix is not a np.ndarray instance.")
        self.density_matrix = density_matrix
    
    def assign_orbital_energies(self, orbital_energies):
        """Assign orbital energies of the wavefunction.

        Parameters
        ----------
        orbital_energies : np.ndarray
            Orbital energies of the wavefunction.

        Raises
        ------
        """
        if orbital_energies is None:
            orbital_energies = np.zeros(self.nspatial)
        elif not isinstance(orbital_energies,np.ndarray):
            raise TypeError("Orbital energies is not a np.ndarray instance.")
        self.orbital_energies = orbital_energies
        
