from moha.system.wavefunction.abc import ABCWaveFunction
import math
import numpy as np

class BaseWaveFunction(ABCWaveFunction):
    """Base wavefunction abstrct class.

    Attributes
    ----------
    nelec : int
        Number of electrons.

    nspatial : int
        Number of spatial orbitals.

    occ : dict
        Occupation number of the wavefunction.

    basis_set
        Basis set of the wavefunction.

    coefficients : np.ndarray
        Coefficientss of the wavefunction.

    Properties
    ----------
    ncoefficients : int
        Number of coefficients.

    nspin : int
        Number of spin orbital.

    spin : int
        Spin of the wavefunction.

    seniority: int
        Seniority of the wavefunction.

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

        basis_set
            Basis set of the wavefunction.

        coefficients : np.ndarray
            Parameters of the wavefunction.
        """
        self.assign_nelec(nelec)
        self.assign_nspatial(nspatial)
        self.assign_occ(occ)
        self.assign_basis_set(basis_set)
        self.assign_coefficients(coefficients)

    @property
    def nspin(self):
        """Return the number of spin orbitals.
       
        Returns
        -------
        nspatial : int
            Number of spin orbitals.

        """
        return self.nspatial*2

    @property
    def ncoefficients(self):
        """Return the number of wavefunction coefficients.
        
        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError        

    @property
    def spin(self):
        """Return the spin of the wavefunction.

        .. math::

            \frac{1}{2}(N_\alpha - N_\beta)
        
        Returns
        -------
        spin : float
            Spin of the wavefunction.

        Notes
        -----
        `None` means that all possible spins are allowed.

        """
        return None

    @property
    def seniority(self):
        """Return the seniority of the wavefunction.

        Seniority of a Slater determinant is its number of unpaired electrons. The seniority of the
        wavefunction is the expected number of unpaired electrons.
        
        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError        

    def assign_nelec(self, nelec):
        """Assign the number of electrons.

        Parameters
        ----------
        nelec : int
            Number of electrons.
            
        Raises
        ------
        TypeError
            If number of electrons is not an integer.
        ValueError
            If number of electrons is not a positive number.

        """
        if not isinstance(nelec, int):
            raise TypeError("Number of electrons must be an integer")
        if nelec <= 0:
            raise ValueError("Number of electrons must be a positive integer")
        self.nelec = nelec

    def assign_nspatial(self, nspatial):
        """Assign the number of spatial orbitals.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbitals.
        
        Raises
        ------
        TypeError
            If number of electrons is not an integer.
        ValueError
            If number of electrons is not a positive number.
        """
        if not isinstance(nspatial, int):
            raise TypeError("Number of spatial orbitals must be an integer")
        if nspatial <= 0:
            raise ValueError("Number of spatial orbitals must be a positive integer")
        self.nspatial = nspatial

    def assign_occ(self, occ):
        """Assign the occupation number of the wavefunction.

        Parameters
        ----------
        occ : dict
            Occupation number of the wavefunction.

        Raises
        ------
        TypeError
            If number of electrons is not an integer.
        ValueError
            If number of electrons is not a positive number.
        """
        if occ == {}:
            occ['alpha'] = math.ceil(self.nelec/2)
            occ['beta'] = math.floor(self.nelec/2)
        elif not isinstance(occ, dict):
            raise TypeError("Occupation number must be dictionary")
        elif not (
            len(occ)==2
            and isinstance(occ['alpha'],int)
            and isinstance(occ['beta'],int)
        ):
            raise ValueError("Occupation number must be a length two dictionary")
        self.occ = occ

    def assign_basis_set(self, basis_set):
        """Assign the basis_set of the wavefunction.

        Parameters
        ----------
        basis_set
            Basis set of the wavefunction.

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError

    def assign_coefficients(self, coefficients):
        """Assign the coefficients of the wavefunction.

        Parameters
        ----------
        coefficients
            Parameters of the wavefunction.
        
        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError
