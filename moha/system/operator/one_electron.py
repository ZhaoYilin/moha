from moha.system.operator.base import BaseOperator
import numpy as np

class RestrictedOneElectronOperator(BaseOperator):
    """Restricted one electron operator class.

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
    dim : int
        Dimension of the operator.

    npsin : int
        Number of spin orbitals.

    spin_orbital_basis_integral(self)
        Return integral matrix in spin orbitals basis.

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

    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        super().__init__(name,nspatial,nelectron,integral)

    @property
    def spin_orbital_basis_integral(self):
        """Return value matrix in spin orbitals basis.
        """
        nspin = self.nspin
        integral = np.zeros((nspin,nspin))
        for i in range(nspin): 
            for j in range(nspin): 
                if i%2==j%2:
                    integral[i,j]  = self.integral[int(i/2),int(j/2)]
        return integral

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
 
    def basis_transformation(self,matrix):
        """Rotate orbitals with a transformation matrix.

        .. math::

            \widetilde{h}_{ab} &= \sum_{ij} C^\dagger_{ai} h_{ij} C_{jb}\\
            \widetilde{g}_{abcd} &= \sum_{ijkl} C^\dagger_{ai} C^\dagger_{bj} g_{ijkl} C_{kc} C_{ld}

        Parameters
        ----------
        matrix : np.ndarray(K, L)
            Transformation matrix.

        Raises
        ------
        TypeError
            If matrix is not a numpy array.
        ValueError
            If shape of matrix does not match up with the shape of the integrals.

        """
        if not isinstance(matrix, np.ndarray):
            raise TypeError("Transformation matrix must be given as a numpy array.")
        if matrix.ndim != 2:
            raise ValueError("Transformation matrix must be two-dimensional.")
        if matrix.shape[0] != self.integral.shape[0]:
            raise ValueError(
                "Shape of the transformation matrix must match with the shape of the " "integrals."
            )
        
        integral = np.einsum("ij,ia->aj", self.integral, matrix)
        integral = np.einsum("aj,jb->ab", integral, matrix)
        return integral

class UnrestrictedOneElectronOperator(BaseOperator):
    """Unrestricted one electron operator class.

    Attributes
    ----------
    name : str
        Name of the operator.

    nspatial : int
        Number of the spatial orbitals.

    nelectron : int
        Number of electrons.

    integral : 2-dictionary of np.ndarray
        Integral value for the operator.

    Property
    --------
    dim : int
        Dimension of the operator.

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
    def __init__(self,name,nspatial,nelectron=1,integral={}):
        super().__init__(name,nspatial,nelectron,integral)

    @property
    def spin_orbital_basis_integral(self):
        """Return value matrix in spin orbitals basis.
        """
        nspin = self.nspin
        integral = np.zeros((nspin,nspin))
        for i in range(nspin): 
            for j in range(nspin): 
                if i%2==j%2==0:
                    integral[i,j]  = self.integral['alpha'][int(i/2),int(j/2)]
                elif i%2==j%2==1:
                    integral[i,j]  = self.integral['beta'][int(i/2),int(j/2)]
        return integral

    def assign_integral(self,integral):
        """Assign integral value for the unrestricted operator.

        Parameters
        ----------
        integral :
            Integral value for the operator.

        Raises
        ------
        TypeError
            If value is not np.ndarray instance.
        """
        if not isinstance(value, dict):
            raise TypeError("Integral must be dict")
        self.integral = integral

    def basis_transformation(self,matrix):
        """Rotate orbitals with a transformation matrix.

        .. math::

            \widetilde{h}_{ab} &= \sum_{ij} C^\dagger_{ai} h_{ij} C_{jb}\\
            \widetilde{g}_{abcd} &= \sum_{ijkl} C^\dagger_{ai} C^\dagger_{bj} g_{ijkl} C_{kc} C_{ld}

        Parameters
        ----------
        matrix : np.ndarray(K, L), 2-dictionary of np.ndarray(K, L)
            Transformation matrix.
            If two transformation matrices are given, the first transforms the alpha orbitals, and
            the second transforms the beta orbitals.
            If only one transformation matrix is given, then the same transformation is applied to
            both the alpha and beta orbitals.        

        Raises
        ------
        TypeError
            If matrix is not a numpy array.
        ValueError
            If shape of matrix does not match up with the shape of the integrals.

        """
        if not isinstance(matrix, np.ndarray):
            raise TypeError("Transformation matrix must be given as a numpy array.")
        if matrix.ndim != 2:
            raise ValueError("Transformation matrix must be two-dimensional.")
        if matrix.shape[0] != self.one_int.shape[0]:
            raise ValueError(
                "Shape of the transformation matrix must match with the shape of the " "integrals."
            )
        for spin in ['alpha','beta']:
            value[spin] = np.einsum("ij,ia->aj", value, matrix[spin])
            value[spin] = np.einsum("aj,jb->ab", value, matrix[spin])
        return value

class RestrictedOverlapOperator(RestrictedOneElectronOperator):
    """Restriced overlap operator class.

    Attributes
    ----------
    name : str
        Name of the operator.

    nspatial : int
        Number of the spatial orbitals.

    nelectron : int
        Number of electrons.

    value : (int,float,np.ndarray)
        Value of the operator.

    Property
    --------
    dim : int
        Dimension of the operator.

    npsin : int
        Number of spin orbitals.

    Methods
    -------
    __init__(self,name,nspatial,value)
        Initialize the Operator.

    build(cls,molecule,basis_set)
        Build the operator instance.

    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.

    assign_nelectron(self,nelectron)
        Assign number of elctrons involved in the integral.

    assign_value(self,value)
        Assign integral value for the operator.
    """
    def __init__(self,name,nspatial,nelectron=1,value=None):
        super().__init__(name,nspatial,nelectron,value)

    @classmethod
    def build(cls,basis_set):
        from moha.system.integral.overlap import overlap
        nspatial = basis_set.size
        value = np.zeros((nspatial,nspatial))     
        for i,oi in enumerate(basis_set.basis):
            for j,oj in enumerate(basis_set.basis):
                if i>=j:
                    value[i][j] = overlap(oi,oj)
        value = value + value.T - np.diag(value.diagonal())
        operator = cls('overlap',nspatial,1,value)
        return operator

