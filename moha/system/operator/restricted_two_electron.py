from moha.system.operator.base import BaseOperator
from moha.system.molecule import Molecule
from moha.system.basis_set import GeneralBasisSet
import numpy as np

class RestrictedTwoElectronOperator(BaseOperator):
    """Restricted two electron operator class.

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
    
    dtype : 
        Data type of the integrals.

    spin_orbital_basis_integral 
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
    def __init__(self,name,nspatial,nelectron=2,integral=None):
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
        """Data type of the integrals.

        Returns
        -------
        dtype : float,int
            Data type of the integrals.

        """
        return self.integral.dtype

    @property
    def spin_orbital_basis_integral(self):
        """Return value matrix in spin orbitals basis.
        """    
        nspin = self.nspin
        integral = np.zeros((nspin,nspin,nspin,nspin))
        for i in range(nspin): 
            for j in range(nspin): 
                for k in range(nspin): 
                    for l in range(nspin): 
                        if i%2==k%2 and j%2==l%2:
                            integral[i,j,k,l]  = self.integral[int(i/2),int(k/2),int(j/2),int(l/2)]
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
        self.integral = np.einsum("ijkl,ia->ajkl", self.integral, matrix)
        self.integral = np.einsum("ajkl,jb->abkl", self.integral, matrix)
        self.integral = np.einsum("abkl,kc->abcl", self.integral, matrix)
        self.integral = np.einsum("abcl,ld->abcd", self.integral, matrix)            

class RestrictedElectronRepulsionOperator(RestrictedTwoElectronOperator):
    """Restricted electron repulsion operator class.

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
    
    basis_transformation(self,matrix)
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self,name,nspatial,nelectron=2,integral=None):
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
        basis_set : GeneralBasisSet
            General basis set instance.

        Raises
        ------
        TypeError
            If basis_set parameter is not a GeneralBasisSet instance.
        """        
        from moha.system.integral.electron_repulsion import electron_repulsion
        from moha.system.auxiliary import eint
        if not isinstance(basis_set,GeneralBasisSet):
            raise TypeError("basis_set parameter must be a GeneralBasisSet instance")       
        nspatial = basis_set.size
        integral_list = []     
        integral = np.zeros((nspatial,nspatial,nspatial,nspatial))
        for i,oi in enumerate(basis_set.bases):
            for j,oj in enumerate(basis_set.bases[0:i+1]):
                for k,ok in enumerate(basis_set.bases[0:i+1]):
                    for l,ol in enumerate(basis_set.bases[0:k+1]):
                        ij = i*(i+1)*0.5+j
                        kl = k*(k+1)*0.5+l
                        if ij>=kl:
                            integral_list.append(electron_repulsion(oi,oj,ok,ol))
        for i in range(nspatial):
            for j in range(nspatial):
                for k in range(nspatial):
                    for l in range(nspatial):
                        integral[i,j,k,l] = integral_list[eint(i,j,k,l)]
 
        operator = cls('electron_repulsion',nspatial,2,integral)
        return operator

    @property
    def coulomb(self):
        """Coulomb integral.

        Returns
        -------
        integral : np.ndarray
            Coulomb integral.        
        """
        nspatial = self.nspatial
        integral = np.zeros((nspatial,nspatial))
        for i in range(nspatial):
            for j in range(nspatial):
                integral[i,j] = self.integral[i,i,j,j]
        return integral

    @property
    def exchange(self):
        """Exchange integral.        

        Returns
        -------
        integral : np.ndarray
            Exchange integral.        
        """
        nspatial = self.nspatial
        integral = np.zeros((nspatial,nspatial))
        for i in range(nspatial):
            for j in range(nspatial):
                integral[i,j] = self.integral[i,j,i,j]
        return integral

    @property
    def double_bar(self):
        """Antisymmetrized two electron integral.
        
        Returns
        -------
        integral : np.ndarray
            Antisymmetrized two electron integral.        
        """
        nspin = self.nspin
        integral = np.zeros((nspin,nspin,nspin,nspin))
        spin_integral = self.spin_orbital_basis_integral
        for i in range(nspin): 
            for j in range(nspin): 
                for k in range(nspin): 
                    for l in range(nspin): 
                        integral[i,j,k,l]  = spin_integral[i,j,k,l] - spin_integral[i,j,l,k]
        return integral


