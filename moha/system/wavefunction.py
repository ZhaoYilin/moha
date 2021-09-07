from abc import ABC,abstractmethod
import numpy as np 

class BaseWaveFunction(ABC):
    """Base wavefunction abstrct class.

    Attributes
    ----------
    nelec : int
        Number of electrons.
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
    assign_basis_set(self, basis_set)
        Assign basis set of the wavefunction.
    assign_coefficients(self, coefficients)
        Assign coefficients of the wavefunction.

    """

    def __init__(self,nelec,nspatial,basis_set=None,coefficients=None):
        """Initialize the wavefunction.

        Parameters
        ----------
        nelec : int
            Number of electrons.
        nspin : int
            Number of spin orbitals.
        dtype : {float, complex, np.float64, np.complex128, None}
            Numpy data type.
            Default is `np.float64`.
        memory : {float, int, str, None}
            Memory available for the wavefunction.
            Default does not limit memory usage (i.e. infinite).

        """
        self.assian_nelec(nelec)
        self.assian_nspatial(nspatial)
        self.assian_basis_set(basis_set)
        self.assian_coefficients(coefficients)

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

        Returns
        -------
        nparams : int
            Number of coefficients.

        """
        return None

    @property
    def spin(self):
        r"""Return the spin of the wavefunction.

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

        Returns
        -------
        seniority : int
            Seniority of the wavefunction.

        Notes
        -----
        `None` means that all possible seniority are allowed.

        """
        return None

    @abstractmethod           
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

    @abstractmethod           
    def assign_nspatial(self, nspatial):
            """Assign the number of spatial orbitals.

            Parameters
            ----------
            nelec : int
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
            if nelec <= 0:
                raise ValueError("Number of spatial orbitals must be a positive integer")
            self.nspatial = nspatial

    @abstractmethod           
    def assign_basis_set(self, basis_set):
        """Assign the basis_set of the wavefunction.

        Parameters
        ----------
        basis_set : Basis class
            basis set of the wavefunction.

        Notes
        -----

        """
        return None

    @abstractmethod           
    def assign_coefficients(self, coefficients):
        """Assign the coefficients of the wavefunction.

        Parameters
        ----------
        params : {np.ndarray, None}
            Parameters of the wavefunction.

        Notes
        -----

        """
        return None
