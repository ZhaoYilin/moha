from moha.system.operator.operator import Operator
from moha.system.molecule import Molecule
from moha.system.basis.basis import Orbital
from moha.system.basis.basis_set import BasisSet

import numpy as np

__all__ = ['NuclearRepulsion','Overlap','Kinetic','NuclearAttraction','ElectronRepulsion']

class NuclearRepulsion(Operator):
    """Nuclear repulsion operator class.
    
    Methods
    -------
    """
    @classmethod
    def build(cls,mol):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        mol : Molecule 
            Molecule instance.

        Raises
        ------
        TypeError
            If parameter mol is not instance of Molecule.

        Return
        ------
        operator: NuclearRepulsion
            Return NuclearRepulsion operator instance.
        """
        if not isinstance(mol, Molecule):
            raise TypeError("Parameter mol must be instance of Molecule.")

        e_nuc = mol.e_nuc
        operator = cls([e_nuc])

        return operator

    def basis_transform(self,matrix):
        """Rotate orbitals with a transformation matrix.

        Parameters
        ----------
        matrix : np.ndarray(K, L)
            Transformation matrix.

        Return
        ------
        NotImplemented
        """
        return NotImplemented
    
class Overlap(Operator):
    """Overlap operator class.
    
    Methods
    -------
    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __call__(self, orb1, orb2):
        """Evaluates overlap between two orbitals.

        Parameters
        ----------
        orb1: Orbital
            The first orbital.

        orb2: Orbital
            The second orbital.

        Return
        ------
        vlaue: float
            Integral value.
        """
        if not isinstance(orb1,Orbital):
            raise TypeError('input must be Orbital instances')

        from moha.system.operator.integral.overlap import Sxyz
        value = 0.0
        for a,cga in enumerate(orb1):
            for b,cgb in enumerate(orb2):
                for i, ci in enumerate(cga.coefs):
                    for j, cj in enumerate(cgb.coefs):
                        value += cga.norm[i]*cgb.norm[j]*ci*cj*\
                            Sxyz(cga.exps[i],cga.shell,cga.origin,
                            cgb.exps[j],cgb.shell,cgb.origin)
        return value

    @classmethod
    def build(cls,basis_set):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : BasisSet
            Hatree Fock basis set instance.

        Raises
        ------
        TypeError
            If basis_set parameter is not a BasisSet instance.

        Return
        ------
        operator: Overlap
            Return an operator instance.
        """
        if not isinstance(basis_set, BasisSet):
            raise TypeError("basis_set parameter must be a BasisSet instance.")
        norb = len(basis_set)    
        operator = np.zeros((norb,norb)).view(cls)
        for i,oi in enumerate(basis_set):
            for j,oj in enumerate(basis_set):
                if i>=j:
                     operator[i,j] = operator(oi,oj)
        operator = operator + operator.T - np.diag(operator.diagonal())
        return operator

class Kinetic(Operator):
    """Kinetic operator class.
    
    Methods
    -------
    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __call__(self, orb1, orb2):
        """Evaluates overlap between two orbitals.

        Parameters
        ----------
        orb1: Orbital
            The first orbital.

        orb2: Orbital
            The second orbital.

        Return
        ------
        vlaue: float
            Integral value.
        """
        if not isinstance(orb1,Orbital):
            raise TypeError('input must be Orbital instances.')

        from moha.system.operator.integral.kinetic import Txyz
        value = 0.0
        for a,cga in enumerate(orb1):
            for b,cgb in enumerate(orb2):
                for i, ci in enumerate(cga.coefs):
                    for j, cj in enumerate(cgb.coefs):
                        value += cga.norm[i]*cgb.norm[j]*ci*cj*\
                            Txyz(cga.exps[i],cga.shell,cga.origin,
                            cgb.exps[j],cgb.shell,cgb.origin)
        return value

    @classmethod
    def build(cls,basis_set):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : BasisSet
            Hatree Fock basis set instance.

        Raises
        ------
        TypeError
            If basis_set parameter is not a BasisSet instance.

        Return
        ------
        operator: Kinetic
            Return an operator instance.
        """
        if not isinstance(basis_set, BasisSet):
            raise TypeError("basis_set parameter must be a BasisSet instance")
        norb = len(basis_set)    
        operator = np.zeros((norb,norb)).view(cls)
        for i,oi in enumerate(basis_set):
            for j,oj in enumerate(basis_set):
                if i>=j:
                     operator[i,j] = operator(oi,oj)
        operator = operator + operator.T - np.diag(operator.diagonal())
        return operator

class NuclearAttraction(Operator):
    """Nuclear attraction operator class.
    
    Methods
    -------
    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __call__(self, orb1, orb2, coordinate):
        """Evaluates overlap between two orbitals.

        Parameters
        ----------
        orb1: Orbital
            The first orbital.

        orb2: Orbital
            The second orbital.
    
        coordinate: 3-list
            Center of nucleus.

        Return
        ------
        vlaue: float
            Integral value.
        """
        if not isinstance(orb1,Orbital):
            raise TypeError('input must be Orbital instances.')

        from moha.system.operator.integral.nuclear_attraction import Vxyz
        value = 0.0
        for a,cga in enumerate(orb1):
            for b,cgb in enumerate(orb2):
                for i, ci in enumerate(cga.coefs):
                    for j, cj in enumerate(cgb.coefs):
                        value += cga.norm[i]*cgb.norm[j]*ci*cj*\
                            Vxyz(cga.exps[i],cga.shell,cga.origin,
                            cgb.exps[j],cgb.shell,cgb.origin,coordinate)
        return value

    @classmethod
    def build(cls,molecule,basis_set):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : BasisSet
            Hatree Fock basis set instance.

        Raises
        ------
        TypeError
            If basis_set parameter is not a BasisSet instance.

        Return
        ------
        operator: Kinetic
            Return an operator instance.
        """
        if not isinstance(molecule, Molecule):
            raise TypeError("molecule parameter must be a Molecule instance")
        if not isinstance(basis_set, BasisSet):
            raise TypeError("basis_set parameter must be a BasisSet instance")
        norb = len(basis_set)    
        operator = np.zeros((norb,norb)).view(cls)
        for i,oi in enumerate(basis_set):
            for j,oj in enumerate(basis_set):
                if i>=j:
                    for atom in molecule:
                        operator[i,j] += -1*atom.number*operator(oi,oj,atom.coordinate)
        operator = operator + operator.T - np.diag(operator.diagonal())
        return operator

class ElectronRepulsion(Operator):
    """Electron repulsion operator class.
    
    Methods
    -------
    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __call__(self, orb1, orb2, orb3, orb4):
        """Evaluates overlap between two orbitals.

        Parameters
        ----------
        orb1: Orbital
            The first orbital.

        orb2: Orbital
            The second orbital.

        orb3: Orbital
            The third orbital.

        orb4: Orbital
            The fourth orbital.

        Return
        ------
        vlaue: float
            Integral value.
        """
        if not isinstance(orb1,Orbital):
            raise TypeError('input must be Orbital instances.')

        from moha.system.operator.integral.electron_repulsion import gxyz
        value = 0.0
        for a,cga in enumerate(orb1):
            for b,cgb in enumerate(orb2):
                for c,cgc in enumerate(orb3):
                    for d,cgd in enumerate(orb4):
                        for i, ci in enumerate(cga.coefs):
                            for j, cj in enumerate(cgb.coefs):
                                for k, ck in enumerate(cgc.coefs):
                                    for l, cl in enumerate(cgd.coefs):
                                        value += cga.norm[i]*cgb.norm[j]*cgc.norm[k]*cgd.norm[l]*\
                                            ci*cj*ck*cl*\
                                            gxyz(cga.exps[i],cga.shell,cga.origin,\
                                            cgb.exps[j],cgb.shell,cgb.origin,\
                                            cgc.exps[k],cgc.shell,cgc.origin,\
                                            cgd.exps[l],cgd.shell,cgd.origin)
        return value

    @classmethod
    def build(cls,basis_set):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        basis_set : BasisSet
            Hatree Fock basis set instance.

        Raises
        ------
        TypeError
            If basis_set parameter is not a BasisSet instance.

        Return
        ------
        operator: Operator
            Return an operator instance.
        """
        from moha.system.operator.integral.auxiliary import eint
        if not isinstance(basis_set, BasisSet):
            raise TypeError("basis_set parameter must be a BasisSet instance")

        norb = len(basis_set)
        integral_list = []    
        operator = np.zeros((norb,norb,norb,norb)).view(cls)
        
        # Calculate integral list with 8-fold permutational symmetry
        for i,oi in enumerate(basis_set):
            for j,oj in enumerate(basis_set[0:i+1]):
                for k,ok in enumerate(basis_set[0:i+1]):
                    for l,ol in enumerate(basis_set[0:k+1]):
                        ij = i*(i+1)*0.5+j
                        kl = k*(k+1)*0.5+l
                        if ij>=kl:
                            integral_list.append(operator(oi,oj,ok,ol))

        # Build integral tensor from list 
        for i in range(norb):
            for j in range(norb):
                for k in range(norb):
                    for l in range(norb):
                        operator[i,j,k,l] = integral_list[eint(i,j,k,l)]

        return operator

    @property
    def coulomb(self):
        """Coulomb integral over spatial orbitals.
        Jij = (ii|jj)

        Returns
        -------
        Jij : np.ndarray
            Coulomb integral.        
        """
        Jij = np.einsum("iijj->ij", self)
        return Jij 

    @property
    def exchange(self):
        """Exchange integral over spatial orbitals.        
        Kij = (ij|ij)

        Returns
        -------
        Kij : np.ndarray
            Exchange integral.        
        """
        Kij = np.einsum("ijij->ij", self)
        return Kij

    def basis_transform(self,C):
        """Transform Eri from atomic orbitals to molecular orbitals.

        Parameters
        ----------
        C : np.ndarray
            Coefficient matrix.

        Raises
        ------
        TypeError
            If coefficient matrix is not a numpy array.
        ValueError
            If shape of coefficient matrix does not match up with the shape of
            the integrals.
        """
        if not isinstance(C, np.ndarray):
            raise TypeError("Coefficient matrix must be given as a numpy array.")
        if C.ndim != 2:
            raise ValueError("Coefficient matrix must be two-dimensional.")
        if C.shape[0] != self.shape[1]:
            raise ValueError(
                "Shape of the coefficient matrix must match with the shape of Eri."
            )
        Eri = np.einsum("ijkl,ia->ajkl", self, C)
        Eri = np.einsum("ajkl,jb->abkl", Eri, C)
        Eri = np.einsum("abkl,kc->abcl", Eri, C)
        Eri = np.einsum("abcl,ld->abcd", Eri, C)
        Eri = Eri.view(self.__class__)
        return Eri
