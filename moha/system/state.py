from abc import ABC
from abc import abstractmethod
import numpy as np 

class WaveFunction(ABC):
    def __init__(self,nelec,nspatial,basis_set=None,coefficients=None):
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
        return np.prod(self.coefficients_shape)

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
           
    def assign_coefficients(self, coefficients):
        """Assign the coefficients of the wavefunction.

        Parameters
        ----------
        params : {np.ndarray, None}
            Parameters of the wavefunction.
        add_noise : bool
            Flag to add noise to the given parameters.

        Raises
        ------
        TypeError
            If `params` is not a numpy array.
            If `params` does not have data type of `float`, `complex`, `np.float64` and
            `np.complex128`.
            If `params` has complex data type and wavefunction has float data type.
        ValueError
            If `params` does not have the same shape as the template_params.

        Notes
        -----
        Depends on dtype, template_params, and nparams.

        """
        if params is None:
            params = self.template_params

        # check if numpy array and if dtype is one of int, float, or complex
        super().assign_params(params)
        params = self.params

        # check shape and dtype
        if params.size != self.nparams:
            raise ValueError("There must be {0} parameters.".format(self.nparams))
        if params.dtype in (complex, np.complex128) and self.dtype in (float, np.float64):
            raise TypeError(
                "If the parameters are `complex`, then the `dtype` of the wavefunction "
                "must be `np.complex128`"
            )
        if params.dtype not in [float, np.float64, complex, np.complex128]:
            raise TypeError("If the parameters are neither float or complex.")

        if len(params.shape) == 1:
            params = params.reshape(self.params_shape)
        elif params.shape != self.params_shape:
            raise ValueError(
                "Parameters must either be flat or have the same shape as the "
                "template, {0}.".format(self.params_shape)
            )

        self.coefficients = coefficients

