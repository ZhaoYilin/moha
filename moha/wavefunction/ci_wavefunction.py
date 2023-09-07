from moha.wavefunction.wavefunction import WaveFunction

import numpy as np 

class CIWaveFunction(WaveFunction):
    """Configuration interaction wavefunction class.

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
        Coefficientss of the wavefunction.

    Properties
    ----------
    ncoefficients : int
        Number of coefficients.
    
    nspin : int
        Number of spin orbital
    
    spin : int
        Spin of the wavefunction
    
    seniority: int
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
    """

    def __init__(self,nelec,nspatial,occ={},basis_set=None,coefficients=None):
        """Initialize the wavefunction.

        Parameters
        ----------
        nelec : int
            Number of electrons.
        
        nspin : int
            Number of spin orbitals.
        
        occ : dict
            Occupation number of the wavefunction.

        dtype : {float, complex, np.float64, np.complex128, None}
            Numpy data type.
            Default is `np.float64`.
        
        memory : {float, int, str, None}
            Memory available for the wavefunction.
            Default does not limit memory usage (i.e. infinite).

        """
        super().__init__(nelec,nspatial,occ,basis_set,coefficients)

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
        return None
    
    def assign_basis_set(self, basis_set):
        """Assign the basis_set of the wavefunction.

        Parameters
        ----------
        basis_set 
            Basis set of the wavefunction.

        Raises
        ------
        TypeError
            If basis set is not a CIBasisSet instance.        
        """
        from moha.posthf.ci.basis_set import CIBasisSet
        if not isinstance(basis_set,CIBasisSet):
            raise TypeError("Basis set must be CIBasisSet instance.")
        self.basis_set = basis_set

    def assign_coefficients(self, coefficients):
        """Assign the coefficients of the wavefunction.

        Parameters
        ----------
        coefficients
            Parameters of the wavefunction.

        Raises
        ------
        TypeError
            If coefficients is not a np.ndarray.        
        """
        if coefficients is None:
            coefficients = np.zeros((self.nspatial))
        elif not isinstance(coefficients,np.ndarray):
            raise TypeError("Coefficients is not a np.ndarray instance.")
        self.coefficients = coefficients
