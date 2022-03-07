from moha.system.wavefunction.base import BaseWaveFunction

import numpy as np 

class HFWaveFunction(BaseWaveFunction):
    """Hatree Fock wavefunction class.

    Attributes
    ----------
    nelec : int
        Number of electrons.

    occ : dict
        Occupation number of the wavefunction.
    
    nspatial : int
        Number of spatial orbitals.
    
    basis_set : Basis
        Basis set of the wavefunction.
    
    coefficients : np.ndarray
        Coefficients of the wavefunction.
    
    density_matrix : np.ndarray
        Density matrix of the wavefunction.
    
    orbital_energies : np.ndarray
        Orbital energies of the wavefunction.

    Properties
    ----------
    ncoefficients : int
        Number of coefficients.
    
    nspin : int
        Number of spin orbital
    
    spin : int
        Spin of the wavefunction
    
    seniority : int
        Seniority of the wavefunction

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

    def __init__(self,nelec,nspatial,occ={},basis_set=None,coefficients=None,density_matrix=None,orbital_energies=None):
        """Initialize the wavefunction.

        Parameters
        ----------
        nelec : int
            Number of electrons.
        
        nspin : int
            Number of spin orbitals.
        
        occ : dict
            Occupation number of the wavefunction.

        basis_set
            Basis set of the wavefunction.

        coefficients : np.ndarray
            Parameters of the wavefunction.
        
        density_matrix : np.ndarray
            Density matrix of the wavefunction.
        
        orbital_energies : np.ndarray
            Orbital energies of the wavefunction.
        """
        super().__init__(nelec,nspatial,occ)
        self.assign_basis_set(basis_set)
        self.assign_coefficients(coefficients)
        self.assign_orbital_energies(orbital_energies)

    @property
    def configuration(self):
        """Return the cofiguration of the wavefunction.

        Returns
        -------
        c : dict
            Configuration of the wavefunction.
        """
        c = {}
        for spin in self.occ:
            c[spin] = [1]*self.occ[spin] + [0]*(self.nspatial - self.occ[spin])
        return c
    
    @property
    def ncoefficients(self):
        """Return the number of wavefunction coefficients.

        Returns
        -------
        ncoefficients : int
            Number of coefficients.
        
        Raises
        ------
        TypeError
            If coefficients is not a np.ndarray instance.
        """
        if not isinstance(self.coefficients,np.ndarray):
            raise TypeError("Coefficients is not a np.ndarray instance.")
        return self.coefficients.size
    
    @property
    def seniority(self):
        """Return the seniority of the wavefunction.

        Seniority of a Slater determinant is its number of unpaired electrons. The seniority of the
        wavefunction is the expected number of unpaired electrons.

        Returns
        -------
        seniority : int
            Seniority of the wavefunction.

        Notes
        -----
        `None` means that all possible seniority are allowed.

        """
        s = sum(map(abs,self.occ['alpha'] - self.occ['beta']))
        return s
    
    def assign_basis_set(self, basis_set):
        """Assign the basis_set of the wavefunction.

        Parameters
        ----------
        basis_set
            Basis set of the wavefunction.

        Raises
        ------
        """
        self.basis_set = basis_set

    def assign_coefficients(self, coefficients):
        """Assign the coefficients of the wavefunction.

        Parameters
        ----------
        coefficients : np.ndarray
            Parameters of the wavefunction.

        Raises
        ------
        """
        if coefficients is None:
            coefficients = np.zeros((self.nspatial,self.nspatial))
        elif not isinstance(coefficients,np.ndarray):
            raise TypeError("Coefficients is not a np.ndarray instance.")
        self.coefficients = coefficients
    
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
        
