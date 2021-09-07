from moha.system.operator.base import BaseOperator
from moha.system.molecule import Molecule
from moha.system.basis import BasisSet
import numpy as np


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

    integral : np.ndarray
        Integral of the operator.

    Property
    --------
    npsin : int
        Number of spin orbitals.

    spin_orbital_basis_integral(self)
        Return integral matrix in spin orbitals basis.

    Methods
    -------
    __init__(self,name,nspatial,integral={})
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
    def __init__(self,name,nspatial,nelectron=1,integral={}):
        """Initialize the Operator.

        Parameters
        ----------
        name : str
            Name of the oprerator.

        nspatial : int
            Number of the spatial orbitals for the operator

        nelectron : int
            Number of elctrons involved in the integral.

        integral : 2-dictionary of np.ndarray
            Integral value for the operator.            
        """
        super().__init__(name,nspatial,nelectron,integral)

    @property
    def spin_orbital_basis_integral(self):
        """Return value matrix in spin orbitals basis.

        Raises
        ------
        ValueError
            If integral not assigned.
        """
        if self.integral == {}:
            raise ValueError("Integral must be assigned")
        nspin = self.nspin
        integral = {'alpha':np.zeros((nspin,nspin)),'beta':np.zeros((nspin,nspin))}
        for i in range(nspin): 
            for j in range(nspin): 
                if i%2==j%2==0:
                    integral['alpha'][i,j]  = self.integral['alpha'][int(i/2),int(j/2)]
                elif i%2==j%2==1:
                    integral['beta'][i,j]  = self.integral['beta'][int(i/2),int(j/2)]
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
        if integral == {}:
            self.integral = integral
        if not (
            isinstance(integral, dict)
            and len(integral) == 2
            and all(isinstance(integral[i], np.ndarray) and integral[i].dtype in [float, complex] for i in integral)
        ):
            raise TypeError(
                "One-electron integrals must be given as a dictionary of two numpy "
                "arrays (with dtype float/complex)."
            )        
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
        if self.integral == {}:
            raise ValueError("Integral must be assigned.")
        if not isinstance(matrix, np.ndarray):
            raise TypeError("Transformation matrix must be given as a numpy array.")
        if matrix.ndim != 2:
            raise ValueError("Transformation matrix must be two-dimensional.")
        if matrix.shape[0] != self.integral.shape[0]:
            raise ValueError(
                "Shape of the transformation matrix must match with the shape of the " "integrals."
            )
        self.integral['alpha'] = np.einsum("ij,ia->aj", self.integral['alpha'], matrix['alpha'])
        self.integral['alpha'] = np.einsum("aj,jb->ab", self.integral['alpha'], matrix['alpha'])

        self.integral['beta'] = np.einsum("ij,ia->aj", self.integral['beta'], matrix['beta'])
        self.integral['beta'] = np.einsum("aj,jb->ab", self.integral['beta'], matrix['beta'])

class UnrestrictedOverlapOperator(UnrestrictedOneElectronOperator):
    """Unrestriced overlap operator class.

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
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set):
        from moha.system.integral.overlap import overlap
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))     
        for i,oi in enumerate(basis_set.basis):
            for j,oj in enumerate(basis_set.basis):
                if i>=j:
                    integral[i][j] = overlap(oi,oj)
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('overlap',nspatial,1,{'alpha':integral,'beta':integral})
        return operator

class UnrestrictedKineticOperator(UnrestrictedOneElectronOperator):
    """Unrestriced kinetic operator class.

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
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set):
        from moha.system.integral.kinetic import kinetic
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))     
        for i,oi in enumerate(basis_set.basis):
            for j,oj in enumerate(basis_set.basis):
                if i>=j:
                    integral[i][j] = kinetic(oi,oj)
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('kinetic',nspatial,1,{'alpha':integral,'beta':integral})
        return operator

class UnrestrictedNuclearAttractionOperator(UnrestrictedOneElectronOperator):
    """Unrestriced nuclear attraction operator class.

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
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,molecule,basis_set):
        from moha.system.integral.nuclear_attraction import nuclear_attraction
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))     
        for i,oi in enumerate(basis_set.basis):
            for j,oj in enumerate(basis_set.basis):
                if i>=j:
                    for atom in molecule.atoms:
                        integral[i][j] += -1*atom.number*nuclear_attraction(oi,oj,atom.coordinate)
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('nuclear_attraction',nspatial,1,{'alpha':integral,'beta':integral})
        return operator

class UnrestrictedMultipoleMomentOperator(UnrestrictedOneElectronOperator):
    """Unrestriced multipole moment operator class.

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
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set,*args):
        from moha.system.integral.multipole_moment import multipole_moment
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))     
        for i,oi in enumerate(basis_set.basis):
            for j,oj in enumerate(basis_set.basis):
                if i>=j:
                    integral[i][j] = multipole_moment(oi,oj,args[0],args[1])
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('multipole_moment',nspatial,1,{'alpha':integral,'beta':integral})
        return operator

class UnrestrictedDifferentialOperator(UnrestrictedOneElectronOperator):
    """Unrestriced multipole moment operator class.

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
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set,*args):
        from moha.system.integral.differential import differential
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))
        for i,oi in enumerate(basis_set.basis):
            for j,oj in enumerate(basis_set.basis):
                if i>=j:
                    integral[i][j] = differential(oi,oj,args[0])
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('differential',nspatial,1,{'alpha':integral,'beta':integral})
        return operator

class UnrestrictedLinearMomentumOperator(UnrestrictedOneElectronOperator):
    """Unrestriced linear momentum operator class.

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
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set):
        from moha.system.integral.linear_momentum import linear_momentum
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))
        integral = integral.astype(np.complex)
        for i,oi in enumerate(basis_set.basis):
            for j,oj in enumerate(basis_set.basis):
                if i>=j:
                    integral[i][j] = linear_momentum(oi,oj)
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('linear_momentum',nspatial,1,{'alpha':integral,'beta':integral})
        return operator

class UnrestrictedAngularMomentumOperator(UnrestrictedOneElectronOperator):
    """Unrestriced angular momentum operator class.

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
    """
    def __init__(self,name,nspatial,nelectron=1,integral=None):
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set):
        from moha.system.integral.angular_momentum import angular_momentum
        nspatial = basis_set.size
        integral = np.zeros((nspatial,nspatial))
        integral = integral.astype(np.complex)
        for i,oi in enumerate(basis_set.basis):
            for j,oj in enumerate(basis_set.basis):
                if i>=j:
                    integral[i][j] = angular_momentum(oi,oj)
        integral = integral + integral.T - np.diag(integral.diagonal())
        operator = cls('angular_momentum',nspatial,1,{'alpha':integral,'beta':integral})
        return operator

