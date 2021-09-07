from moha.system.operator.base import BaseOperator
from moha.system.molecule import Molecule
from moha.system.basis import BasisSet
import numpy as np

class UnrestrictedTwoElectronOperator(BaseOperator):
    """Unrestricted two electron operator class.

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

    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self,name,nspatial,nelectron=2,integral=None):
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
        integral = {'alpha':np.zeros((nspin,nspin,nspin,nspin)),'beta':np.zeros((nspin,nspin,nspin,npsin))}
        for i in range(nspin): 
            for j in range(nspin): 
                for k in range(nspin): 
                    for l in range(nspin): 
                        if i%2==k%2==0 and j%2==l%2==0:
                            integral[i,j,k,l]  = self.integral['alpha'][int(i/2),int(k/2),int(j/2),int(l/2)]
                        elif i%2==k%2==1 and j%2==l%2==1:
                            integral[i,j,k,l]  = self.integral['beta'][int(i/2),int(k/2),int(j/2),int(l/2)]
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
        self.integral['alpha'] = np.einsum("ijkl,ia->ajkl", self.integral['alpha'], matrix['alpha'])
        self.integral['alpha'] = np.einsum("ajkl,jb->abkl", self.integral['alpha'], matrix['alpha'])
        self.integral['alpha'] = np.einsum("abkl,kc->abcl", self.integral['alpha'], matrix['alpha'])
        self.integral['alpha'] = np.einsum("abcl,ld->abcd", self.integral['alpha'], matrix['alpha'])            

        self.integral['beta'] = np.einsum("ijkl,ia->ajkl", self.integral['beta'], matrix['beta'])
        self.integral['beta'] = np.einsum("ajkl,jb->abkl", self.integral['beta'], matrix['beta'])
        self.integral['beta'] = np.einsum("abkl,kc->abcl", self.integral['beta'], matrix['beta'])
        self.integral['beta'] = np.einsum("abcl,ld->abcd", self.integral['beta'], matrix['beta']) 

class UnrestrictedElectronRepulsionOperator(UnrestrictedTwoElectronOperator):
    """Unrestricted electron repulsion operator class.

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
    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __init__(self,name,nspatial,nelectron=2,integral=None):
        super().__init__(name,nspatial,nelectron,integral)

    @classmethod
    def build(cls,basis_set):
        from moha.system.integral.electron_repulsion import electron_repulsion
        from moha.system.auxiliary import eint
        nspatial = basis_set.size
        integral_list = []     
        integral = np.zeros((nspatial,nspatial,nspatial,nspatial))
        for i,oi in enumerate(basis_set.basis):
            for j,oj in enumerate(basis_set.basis[0:i+1]):
                for k,ok in enumerate(basis_set.basis[0:i+1]):
                    for l,ol in enumerate(basis_set.basis[0:k+1]):
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
        """
        """
        nspatial = self.nspatial
        integral = np.zeros((nspatial,nspatial))
        for i in range(nspatial):
            for j in range(nspatial):
                integral[i,j] = self.integral[i,i,j,j]
        return integral

    @property
    def exchange(self):
        nspatial = self.nspatial
        integral = np.zeros((nspatial,nspatial))
        for i in range(nspatial):
            for j in range(nspatial):
                integral[i,j] = self.integral[i,j,i,j]
        return integral

    @property
    def double_bar(self):
        """
        antisymmetrized two electron integral
        """
        nspin = self.nspin
        integral = np.zeros((nspin,nspin,nspin,nspin))
        spin_integral = self.spin_orbital_basis_integral
        for i in range(nspin): 
            for j in range(nspin): 
                for k in range(nspin): 
                    for l in range(nspin): 
                        integral[i,j,k,l]  = spin_integral[i,j,k,l] - spin_integral[i,j,l,k]
        return value_prime


