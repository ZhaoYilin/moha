from abc import ABC,abstractmenthod

class BaseHamiltonian(ABC):
    """Hamiltonian for a Schrodinger equation.

    Attributes
    ----------
    nspatial : int
        Number of spatial orbitals.

    operators : dict
        Dictionary of operators for the Hamiltonian

    Methods
    -------
    __init__(self, nspatial)
        Initialize the Hamiltonian.

    assign_operators(self, operator)
        Assigns operator integral to the Hamiltonian.

    """

    def __init__(self, nspatial):
        """Initialize the Hamiltonian.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbitals.

        """
        self.assign_nspatial(nspatial)

    def assign_nspatial(self, nspatial):
        """Assign the number of spatial orbitals.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbital.

        """
        raise NotImplementedError

    @abc.abstractproperty
    def nspin(self):
        """Return the number of spin orbitals.

        Returns
        -------
        nspin : int
            Number of spin orbitals.

        """
        raise NotImplementedError

    def assign_oprators(self, operator):
        """Assign operator.

        Parameters
        ----------
        one_int : np.ndarray(K, K)
            One electron integrals.
        two_int : np.ndarray(K, K, K, K)
            Two electron integrals.
            Uses physicist's notation.

        Raises
        ------
        NotImplementedError
            If called.

        """
        raise NotImplementedError

