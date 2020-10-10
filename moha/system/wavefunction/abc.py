from abc import ABC,abstractmethod,abstractproperty

class ABCWaveFunction(ABC):
    """Abstract wavefunction abstrct class.

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

    Abstract Properties
    -------------------
    ncoefficients : int
        Number of coefficients.

    nspin : int
        Number of spin orbital.

    spin : int
        Spin of the wavefunction.

    seniority: int
        Seniority of the wavefunction.

    Abstract Methods
    ----------------
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
        self.assign_nelec(nelec)
        self.assign_nspatial(nspatial)
        self.assign_occ(occ)
        self.assign_basis_set(basis_set)
        self.assign_coefficients(coefficients)

    @abstractproperty
    def nspin(self):
        """Return the number of spin orbitals.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError

    @abstractproperty
    def ncoefficients(self):
        """Return the number of wavefunction coefficients.

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError

    @abstractproperty
    def spin(self):
        """Return the spin of the wavefunction.

        .. math::

            \frac{1}{2}(N_\alpha - N_\beta)

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError

    @abstractproperty
    def seniority(self):
        """Return the seniority of the wavefunction.

        Seniority of a Slater determinant is its number of unpaired electrons. The seniority of the
        wavefunction is the expected number of unpaired electrons.

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError

    @abstractmethod           
    def assign_nelec(self, nelec):
        """Assign the number of electrons.

        Parameters
        ----------
        nelec : int
            Number of electrons.

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError

    @abstractmethod           
    def assign_nspatial(self, nspatial):
        """Assign the number of spatial orbitals.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbitals.

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError

    @abstractmethod           
    def assign_occ(self, occ):
        """Assign the occupation number of the wavefunction.

        Parameters
        ----------
        occ : dict
            Occupation number of the wavefunction.

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError

    @abstractmethod           
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

    @abstractmethod           
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
