from moha.system.operator.base import OperatorNames, BaseOperator
from moha.system.molecule import Molecule
from moha.system.basis_set.hf_basis_set import HFBasisSet
import numpy as np

class OneElectronOperator(BaseOperator):
    """One electron operator class.

    Attributes
    ----------
    integral : np.ndarray
        Integral of the operator.
    
    name : OperatorNames
        Name of the operator.

    Property
    --------
    npatial : int
        Number of spatial orbitals.

    npsin : int
        Number of spin orbitals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.

    Methods
    -------
    __init__(self, integral, name)
        Initialize the operator

    assign_name(self,name)
        Assign name to the operator.

    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self, integral, name):
        """Initialize the operator

        Parameters
        ----------
        integral : ndarray
            Integral value for the operator.
        
        name : OperatorNames
            Name of the operator.
        """
        super().__init__(integral, name)
    
    @property        
    def nspatial(self):
        """Number of spatial orbitals.
        
        Returns
        -------
        nspatial : int
            Number of spatial orbitals.
        """
        return self.shape[0]
    
    @property        
    def nspin(self):
        """Number of spin orbitals.
        
        Returns
        -------
        nspin : int
            Number of spin orbitals.
        """
        return self.shape[0]*2
    
    @property
    def spin_orbital_basis_integral(self):
        """Return integral in spin orbitals basis.
        
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
                    integral[i,j]  = self[int(i/2),int(j/2)]
        return integral
    
    @classmethod
    def zeros(cls, shape, dtype=float):
        """generate instance with zeros elements ndarray.

        Parameters
        ----------
        shape : tuple/list
            Shape of the ndarray.

        Raises
        ------
        ValueError
            If shape is not two dimension.
        """
        if not len(shape)==2:
            raise ValueError("Shape of the ndarray must be two dimension")
        
        super().zeros(shape, dtype=float)

    @classmethod
    def ones(cls, shape, dtype=float):
        """generate instance with ones elements ndarray.

        Parameters
        ----------
        shape : tuple/list
            Shape of the ndarray.
        
        Raises
        ------
        ValueError
            If shape is not two dimension.
        """
        if not len(shape)==2:
            raise ValueError("Shape of the ndarray must be two dimension")

        super().ones(shape, dtype=float)
    
    @classmethod
    def random(cls, shape, dtype=float):
        """generate instance with random elements ndarray.

        Parameters
        ----------
        shape : tuple/list
            Shape of the ndarray.
        
        Raises
        ------
        ValueError
            If shape is not two dimension.
        """
        if not len(shape)==2:
            raise ValueError("Shape of the ndarray must be two dimension")

        super().random(shape, dtype=float)

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
        if matrix.shape[0] != self.shape[0]:
            raise ValueError(
                "Shape of the transformation matrix must match with the shape of the integrals."
            )
        
        integral = np.einsum("ij,ia->aj", self, matrix)
        integral = np.einsum("aj,jb->ab", integral, matrix)
        self[:] = integral[:]
        return self

class OverlapOperator(OneElectronOperator):
    """Overlap operator class.

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
    __init__(self, integral)
        Initialize the operator

    assign_name(self,name)
        Assign name to the operator.
    
    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self, integral, name=OperatorNames.S):
        """Initialize the operator

        Parameters
        ----------
        integral : ndarray
            Integral value for the operator.
        
        name : OperatorNames
            name of the operator.
        """
        super().__init__(integral, name)

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
        operator = cls(integral,OperatorNames.S)
        return operator

class KineticOperator(OneElectronOperator):
    """Kinetic operator class.

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
    def __init__(self, integral, name=OperatorNames.T):
        """Initialize the operator

        Parameters
        ----------
        integral : ndarray
            Integral value for the operator.
        
        name : OperatorNames
            name of the operator.
        """
        super().__init__(integral, name)

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
        operator = cls(integral,OperatorNames.T)
        return operator

class NuclearAttractionOperator(OneElectronOperator):
    """Nuclear attraction operator class.

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
    def __init__(self, integral, name=OperatorNames.V):
        """Initialize the operator

        Parameters
        ----------
        integral : ndarray
            Integral value for the operator.
        
        name : OperatorNames
            name of the operator.
        """
        super().__init__(integral, name)
    
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
        operator = cls(integral,OperatorNames.V)
        return operator

class MultipoleMomentOperator(OneElectronOperator):
    """Multipole moment operator class.

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
    def __init__(self, integral, name=OperatorNames.MM):
        """Initialize the operator

        Parameters
        ----------
        integral : ndarray
            Integral value for the operator.
        
        name : OperatorNames
            name of the operator.
        """
        super().__init__(integral, name)

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
        operator = cls(integral, OperatorNames.MM)
        return operator

class DifferentialOperator(OneElectronOperator):
    """Differential operator class.

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
    def __init__(self, integral, name=OperatorNames.Diff):
        """Initialize the operator

        Parameters
        ----------
        integral : ndarray
            Integral value for the operator.
        
        name : OperatorNames
            name of the operator.
        """
        super().__init__(integral, name)
    
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
        operator = cls(integral,OperatorNames.Diff)
        return operator

class LinearMomentumOperator(OneElectronOperator):
    """Linear momentum operator class.

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
    def __init__(self, integral, name=OperatorNames.LM):
        """Initialize the operator

        Parameters
        ----------
        integral : ndarray
            Integral value for the operator.
        
        name : OperatorNames
            name of the operator.
        """
        super().__init__(integral, name)
    
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
        operator = cls(integral,OperatorNames.LM)
        return operator

class AngularMomentumOperator(OneElectronOperator):
    """Angular momentum operator class.

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
    def __init__(self, integral, name=OperatorNames.AM):
        """Initialize the operator

        Parameters
        ----------
        integral : ndarray
            Integral value for the operator.
        
        name : OperatorNames
            name of the operator.
        """
        super().__init__(integral, name)
    
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
        operator = cls(integral,OperatorNames.AM)
        return operator
