from moha.system.operator.operator import Operator
from moha.system.basis.basis import Orbital
from moha.system.basis.basis_set import BasisSet

import numpy as np

__all__ = []

class MultipoleMoment(Operator):
    """Multipole moment operator class.

    Property
    --------
    nspatial : int
        Number of the spatial orbitals.

    npsin : int
        Number of spin orbitals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.
    
    Methods
    -------
    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __call__(self, orb1, orb2, lmn, coordinate):
        """Evaluates overlap between two orbitals.

        Parameters
        ----------
        orb1: Orbital
            The first orbital.

        orb2: Orbital
            The second orbital.
        
        lmn: 3-list
            The second orbital.

        coordinate: 3-list
            Coordinate.
            
        Return
        ------
        vlaue: float
            Integral value.
        """
        if not isinstance(orb1,Orbital):
            raise TypeError('input must be Orbital instances.')

        from moha.system.operator.integral.multipole_moment import MMxyz
        value = 0.0
        for a,cga in enumerate(orb1):
            for b,cgb in enumerate(orb2):
                for i, ci in enumerate(cga.coefs):
                    for j, cj in enumerate(cgb.coefs):
                        value += cga.norm[i]*cgb.norm[j]*ci*cj*\
                            MMxyz(cga.exps[i],cga.shell,cga.origin,
                            cgb.exps[j],cgb.shell,cgb.origin,lmn,coordinate)
        return value

    @classmethod
    def build(cls,basis_set, lmn, coordinate):
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
        nspatial = len(basis_set)    
        operator = np.zeros((nspatial,nspatial)).view(cls)
        for i,oi in enumerate(basis_set):
            for j,oj in enumerate(basis_set):
                if i>=j:
                     operator[i,j] = operator(oi,oj,lmn,coordinate)
        operator = operator + operator.T - np.diag(operator.diagonal())
        return operator

class Differential(Operator):
    """Differential operator class.

    Property
    --------
    nspatial : int
        Number of the spatial orbitals.

    npsin : int
        Number of spin orbitals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.
    
    Methods
    -------
    basis_transformation(self,matrix):
        Rotate orbitals with a transformation matrix.
    """
    def __call__(self, orb1, orb2, es):
        """Evaluates overlap between two orbitals.

        Parameters
        ----------
        orb1: Orbital
            The first orbital.

        orb2: Orbital
            The second orbital.
        
        es: 3-list
            The second orbital.
            
        Return
        ------
        vlaue: float
            Integral value.
        """
        if not isinstance(orb1,Orbital):
            raise TypeError('input must be Orbital instances.')

        from moha.system.integral.differential import Dxyz
        value = 0.0
        for a,cga in enumerate(orb1):
            for b,cgb in enumerate(orb2):
                for i, ci in enumerate(cga.coefs):
                    for j, cj in enumerate(cgb.coefs):
                        value += cga.norm[i]*cgb.norm[j]*ci*cj*\
                            Dxyz(cga.exps[i],cga.shell,cga.origin,
                            cgb.exps[j],cgb.shell,cgb.origin,es)
        return value

    @classmethod
    def build(cls,basis_set,es):
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
        nspatial = len(basis_set)    
        operator = np.zeros((nspatial,nspatial)).view(cls)
        for i,oi in enumerate(basis_set):
            for j,oj in enumerate(basis_set):
                if i>=j:
                     operator[i,j] = operator(oi,oj,es)
        operator = operator + operator.T - np.diag(operator.diagonal())
        return operator

class LinearMomentum(Operator):
    """Linear momentum operator class.

    Property
    --------
    nspatial : int
        Number of the spatial orbitals.

    npsin : int
        Number of spin orbitals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.
    
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

        from moha.system.integral.linear_momentum import LMxyz
        value = 0.0
        for a,cga in enumerate(orb1):
            for b,cgb in enumerate(orb2):
                for i, ci in enumerate(cga.coefs):
                    for j, cj in enumerate(cgb.coefs):
                        value += cga.norm[i]*cgb.norm[j]*ci*cj*\
                            LMxyz(cga.exps[i],cga.shell,cga.origin,
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
        nspatial = len(basis_set)    
        operator = np.zeros((nspatial,nspatial)).view(cls)
        for i,oi in enumerate(basis_set):
            for j,oj in enumerate(basis_set):
                if i>=j:
                     operator[i,j] = operator(oi,oj)
        operator = operator + operator.T - np.diag(operator.diagonal())
        return operator

class AngularMomentum(Operator):
    """Angular momentum operator class.

    Property
    --------
    nspatial : int
        Number of the spatial orbitals.

    npsin : int
        Number of spin orbitals.

    spin_orbital_basis_integral : np.ndarray
        Return integral matrix in spin orbitals basis.
    
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

        from moha.system.integral.angular_momentum import AMxyz
        value = 0.0
        for a,cga in enumerate(orb1):
            for b,cgb in enumerate(orb2):
                for i, ci in enumerate(cga.coefs):
                    for j, cj in enumerate(cgb.coefs):
                        value += cga.norm[i]*cgb.norm[j]*ci*cj*\
                            AMxyz(cga.exps[i],cga.shell,cga.origin,
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
        nspatial = len(basis_set)    
        operator = np.zeros((nspatial,nspatial)).view(cls)
        for i,oi in enumerate(basis_set):
            for j,oj in enumerate(basis_set):
                if i>=j:
                     operator[i,j] = operator(oi,oj)
        operator = operator + operator.T - np.diag(operator.diagonal())
        return operator
