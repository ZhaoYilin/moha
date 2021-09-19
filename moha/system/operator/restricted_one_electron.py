from moha.system.operator.base import BaseOperator
from moha.system.molecule import Molecule
from moha.system.basis_set.hf_basis_set import HFBasisSet
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
    npsin : int
        Number of spin orbitals.

    dtype : float64, complex128
        Data type of the integrals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.

    Methods
    -------
    __init__(self,name,nspatial,nelectron=1,integral=None)
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
        """Initialize the operator

        Parameters
        ----------
        name : str
            Name of the operator.

        nspatial : int
            Number of the spatial orbitals for the operator

        nelectron : int
            Number of elctrons involved in the integral.

        integral
            Integral value for the operator.

        """
        super().__init__(name,nspatial,nelectron,integral)

    @property
    def dtype(self):
        """Return the data type of the integrals.

        Returns
        -------
        dtype : float64, complex128
            Data type of the integrals.

        """
        return self.integral.dtype

    @property
    def spin_orbital_basis_integral(self):
        """Return value matrix in spin orbitals basis.
        
        Returns
        -------
        integral : np.ndarray
            Integral of the operator in spin orbital basis.
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
            If integrals are not provided as a numpy array of dtype float/complex.
        ValueError
            If integrals are not given as a (two-dimensional) square matrix.
        """
        if integral is None:
            self.integral = np.zeros((self.nspatial,self.nspatial))
        if not (
            isinstance(integral, np.ndarray)
            and integral.dtype in [float, complex]
        ):
            raise TypeError(
                "Integrals must be given as a numpy array with dtype of float or complex."
            )
        if not (integral.ndim == 2 and integral.shape[0] == integral.shape[1]):
            raise ValueError("One-electron integrals be a (two-dimensional) square matrix.")            
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
                "Shape of the transformation matrix must match with the shape of the integrals."
            )
        
        self.integral = np.einsum("ij,ia->aj", self.integral, matrix)
        self.integral = np.einsum("aj,jb->ab", self.integral, matrix)

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

    integral : np.ndarray
        Integral of the operator.

    Property
    --------
    npsin : int
        Number of spin orbitals.

    dtype : float64, complex128
        Data type of the integrals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.
    
    Methods
    -------
    __init__(self,name,nspatial,nelectron=1,integra=None)
        Initialize the Operator.

    build(cls,molecule,basis_set)
        Build the operator instance.

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
        """Initialize the operator

        Parameters
        ----------
        name : str
            Name of the operator.

        nspatial : int
            Number of the spatial orbitals for the operator

        nelectron : int
            Number of elctrons involved in the integral.

        integral
            Integral value for the operator.
        """
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : HFBasisSet
            Hatree Fock basis set instance.

        Raises
        ------
        TypeError
            If basis_set parameter is not a HFBasisSet instance.
        """
        from moha.system.integral.overlap import overlap
        if not isinstance(basis_set,HFBasisSet):
            raise TypeError("basis_set parameter must be a HFBasisSet instance")    
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))     
        for i,oi in enumerate(basis_set.bases):
            for j,oj in enumerate(basis_set.bases):
                if i>=j:
                    integral[i][j] = overlap(oi,oj)
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('overlap',nspatial,1,integral)
        return operator

class RestrictedKineticOperator(RestrictedOneElectronOperator):
    """Restriced kinetic operator class.

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

    dtype : float64, complex128
        Data type of the integrals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.
    
    Methods
    -------
    __init__(self,name,nspatial,nelectron=1,integra=None)
        Initialize the Operator.

    build(cls,molecule,basis_set)
        Build the operator instance.

    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.

    assign_nelectron(self,nelectron)
        Assign number of elctrons involved in the integral.

    assign_integral(self,integral)
        Assign integral value for the operator.
    
    basis_transformation(self,matrix)
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        """Initialize the operator

        Parameters
        ----------
        name : str
            Name of the operator.

        nspatial : int
            Number of the spatial orbitals for the operator

        nelectron : int
            Number of elctrons involved in the integral.

        integral
            Integral value for the operator.
        """
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : HFBasisSet
            Hatree Fock basis set instance.

        Raises
        ------
        TypeError
            If basis_set parameter is not a HFBasisSet instance.
        """
        from moha.system.integral.kinetic import kinetic
        if not isinstance(basis_set,HFBasisSet):
            raise TypeError("basis_set parameter must be a HFBasisSet instance")    
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))     
        for i,oi in enumerate(basis_set.bases):
            for j,oj in enumerate(basis_set.bases):
                if i>=j:
                    integral[i][j] = kinetic(oi,oj)
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('kinetic',nspatial,1,integral)
        return operator

class RestrictedNuclearAttractionOperator(RestrictedOneElectronOperator):
    """Restriced nuclear attraction operator class.

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

    dtype : float64, complex128
        Data type of the integrals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.
    
    Methods
    -------
    __init__(self,name,nspatial,nelectron=1,integra=None)
        Initialize the Operator.

    build(cls,molecule,basis_set)
        Build the operator instance.

    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.

    assign_nelectron(self,nelectron)
        Assign number of elctrons involved in the integral.

    assign_integral(self,integral)
        Assign integral value for the operator.
    
    basis_transformation(self,matrix)
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        """Initialize the operator

        Parameters
        ----------
        name : str
            Name of the operator.

        nspatial : int
            Number of the spatial orbitals for the operator

        nelectron : int
            Number of elctrons involved in the integral.

        integral
            Integral value for the operator.
        """
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,molecule,basis_set):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        moleucle : Molecule
            Molecule instance.

        basis_set : HFBasisSet
            Hatree Fock basis set instance.

        Raises
        ------
        TypeError
            If molecule parameter is not a Molecule instance.
            If basis_set parameter is not a HFBasisSet instance.
        """
        from moha.system.integral.nuclear_attraction import nuclear_attraction
        if not isinstance(molecule,Molecule):
            raise TypeError("molecule parameter must be a Molecule instance")
        if not isinstance(basis_set,HFBasisSet):
            raise TypeError("basis_set parameter must be a HFBasisSet instance")    
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))     
        for i,oi in enumerate(basis_set.bases):
            for j,oj in enumerate(basis_set.bases):
                if i>=j:
                    for atom in molecule.atoms:
                        integral[i][j] += -1*atom.number*nuclear_attraction(oi,oj,atom.coordinate)
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('nuclear_attraction',nspatial,1,integral)
        return operator

class RestrictedMultipoleMomentOperator(RestrictedOneElectronOperator):
    """Restriced multipole moment operator class.

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

    dtype : float64, complex128
        Data type of the integrals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.
    
    Methods
    -------
    __init__(self,name,nspatial,nelectron=1,integra=None)
        Initialize the Operator.

    build(cls,basis_set,*args)
        Build the operator instance.

    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.

    assign_nelectron(self,nelectron)
        Assign number of elctrons involved in the integral.

    assign_integral(self,integral)
        Assign integral value for the operator.
    
    basis_transformation(self,matrix)
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        """Initialize the operator

        Parameters
        ----------
        name : str
            Name of the operator.

        nspatial : int
            Number of the spatial orbitals for the operator

        nelectron : int
            Number of elctrons involved in the integral.

        integral
            Integral value for the operator.
        """
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set,*args):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : HFBasisSet
            Hatree Fock basis set instance.

        *args
            Arguments for multipole moment

        Raises
        ------
        TypeError
            If basis_set parameter is not a HFBasisSet instance.
        """
        from moha.system.integral.multipole_moment import multipole_moment
        if not isinstance(basis_set,HFBasisSet):
            raise TypeError("basis_set parameter must be a HFBasisSet instance")    
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))     
        for i,oi in enumerate(basis_set.bases):
            for j,oj in enumerate(basis_set.bases):
                if i>=j:
                    integral[i][j] = multipole_moment(oi,oj,args[0],args[1])
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('multipole_moment',nspatial,1,integral)
        return operator

class RestrictedDifferentialOperator(RestrictedOneElectronOperator):
    """Restriced multipole moment operator class.

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

    dtype : float64, complex128
        Data type of the integrals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.
    
    Methods
    -------
    __init__(self,name,nspatial,nelectron=1,integra=None)
        Initialize the Operator.

    build(cls,basis_set,*args)
        Build the operator instance.

    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.

    assign_nelectron(self,nelectron)
        Assign number of elctrons involved in the integral.

    assign_integral(self,integral)
        Assign integral value for the operator.
    
    basis_transformation(self,matrix)
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        """Initialize the operator

        Parameters
        ----------
        name : str
            Name of the operator.

        nspatial : int
            Number of the spatial orbitals for the operator

        nelectron : int
            Number of elctrons involved in the integral.

        integral
            Integral value for the operator.
        """
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set,*args):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : HFBasisSet
            Hatree Fock basis set instance.

        *args
            Arguments for multipole moment.

        Raises
        ------
        TypeError
            If basis_set parameter is not a HFBasisSet instance.
        """
        from moha.system.integral.differential import differential
        if not isinstance(basis_set,HFBasisSet):
            raise TypeError("basis_set parameter must be a HFBasisSet instance")    
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))
        for i,oi in enumerate(basis_set.bases):
            for j,oj in enumerate(basis_set.bases):
                if i>=j:
                    integral[i][j] = differential(oi,oj,args[0])
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('differential',nspatial,1,integral)
        return operator

class RestrictedLinearMomentumOperator(RestrictedOneElectronOperator):
    """Restriced linear momentum operator class.

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

    dtype : float64, complex128
        Data type of the integrals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.
    
    Methods
    -------
    __init__(self,name,nspatial,nelectron=1,integra=None)
        Initialize the Operator.

    build(cls,basis_set,*args)
        Build the operator instance.

    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.

    assign_nelectron(self,nelectron)
        Assign number of elctrons involved in the integral.

    assign_integral(self,integral)
        Assign integral value for the operator.
    
    basis_transformation(self,matrix)
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        """Initialize the operator

        Parameters
        ----------
        name : str
            Name of the operator.

        nspatial : int
            Number of the spatial orbitals for the operator

        nelectron : int
            Number of elctrons involved in the integral.

        integral
            Integral value for the operator.
        """
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : HFBasisSet
            Hatree Fock basis set instance.

        Raises
        ------
        TypeError
            If basis_set parameter is not a HFBasisSet instance.
        """
        from moha.system.integral.linear_momentum import linear_momentum
        if not isinstance(basis_set,HFBasisSet):
            raise TypeError("basis_set parameter must be a HFBasisSet instance")    
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))
        integral = integral.astype(np.complex)
        for i,oi in enumerate(basis_set.bases):
            for j,oj in enumerate(basis_set.bases):
                if i>=j:
                    integral[i][j] = linear_momentum(oi,oj)
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('linear_momentum',nspatial,1,integral)
        return operator

class RestrictedAngularMomentumOperator(RestrictedOneElectronOperator):
    """Restriced angular momentum operator class.

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

    dtype : float64, complex128
        Data type of the integrals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.
    
    Methods
    -------
    __init__(self,name,nspatial,nelectron=1,integra=None)
        Initialize the Operator.

    build(cls,basis_set,*args)
        Build the operator instance.

    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.

    assign_nelectron(self,nelectron)
        Assign number of elctrons involved in the integral.

    assign_integral(self,integral)
        Assign integral value for the operator.
    
    basis_transformation(self,matrix)
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        """Initialize the operator

        Parameters
        ----------
        name : str
            Name of the operator.

        nspatial : int
            Number of the spatial orbitals for the operator

        nelectron : int
            Number of elctrons involved in the integral.

        integral
            Integral value for the operator.
        """
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : HFBasisSet
            Hatree Fock basis set instance.

        Raises
        ------
        TypeError
            If basis_set parameter is not a HFBasisSet instance.
        """
        from moha.system.integral.angular_momentum import angular_momentum
        if not isinstance(basis_set,HFBasisSet):
            raise TypeError("basis_set parameter must be a HFBasisSet instance")    
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))
        integral = integral.astype(np.complex)
        for i,oi in enumerate(basis_set.bases):
            for j,oj in enumerate(basis_set.bases):
                if i>=j:
                    integral[i][j] = angular_momentum(oi,oj)
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('angular_momentum',nspatial,1,integral)
        return operator
