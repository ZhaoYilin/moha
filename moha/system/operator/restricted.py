from moha.system.operator.base import BaseOperator
import numpy as np


class RestrictedOperator(BaseOperator):
    """Restricted operator class.

    Attributes
    ----------
    name : str
        Name of the operator.

    nspatial : int
        Number of the spatial orbitals.

    nelectron : int
        Number of electrons.

    integral : np.ndarray
        Integral of the operator.

    Property
    --------
    npsin : int
        Number of spin orbitals.

    Methods
    -------
    __init__(self,name,nspatial,integral)
        Initialize the Operator.

    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.

    assign_nelectron(self,nelectron)
        Assign number of elctrons involved in the integral.

    assign_integral(self,integral)
        Assign integral value for the operator.
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        super().__init__(name,nspatial,nelectron,integral)

    def assign_integral(self,integral):
        """Assign integral value for the restricted operator.

        Parameters
        ----------
        integral :
            Integral value for the operator.

        Raises
        ------
        TypeError
            If integral is not np.ndarray instance.
        """
        if integral is None:
            self.integral = integral
        elif not isinstance(integral, np.ndarray):
            raise TypeError("Integral must be np.ndarray")
        else:
            self.integral = integral
 
