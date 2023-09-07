import numpy as np
import copy
from moha.basis.gauss import PrimitiveGaussian,ContractedGaussian

__all__ = ['Orbital']

class Orbital(list):
    """Predetermined linear combination of radial parts of GTOs
    Atomic orbtial represented by Contracted Gaussian functions.
    
    Attributes
    ----------
    
    Properties
    ----------
    
    Methods
    -------
    """
    def __call__(self,x,y,z):
        orb = np.zeros(x.shape)
        for cg in self:
            orb += cg(x,y,z) 
        return orb

    def __add__(self,other):
        if isinstance(other,self.__class__):
            orbitals = self.__class__()
            orbitals.extend(copy.deepcopy(self))
            orbitals.extend(copy.deepcopy(other))
            for i,orbi in enumerate(orbitals.copy()):
                for j,orbj in enumerate(orbitals[i+1:].copy()):
                    if orbi.parallel(orbj):
                        coefs = [a+b for a,b in zip(orbi.coefs,orbj.coefs)]
                        orbitals[i].coefs = coefs
                        orbitals.remove(orbj)
            for orb in orbitals.copy():
                if all(orb.coefs)==False:
                    orbitals.remove(orb)
            return orbitals
        else:
            return NotImplemented

    def __sub__(self,other):
        if isinstance(other,self.__class__):
            orbitals = self.__class__()
            orbitals.extend(copy.deepcopy(self))
            orbitals.extend(copy.deepcopy(-other))
            for i,orbi in enumerate(orbitals.copy()):
                for j,orbj in enumerate(orbitals[i+1:].copy()):
                    if orbi.parallel(orbj):
                        coefs = [a+b for a,b in zip(orbi.coefs,orbj.coefs)]
                        orbitals[i].coefs = coefs
                        orbitals.remove(orbj)
            for orb in orbitals.copy():
                if all(orb.coefs)==False:
                    orbitals.remove(orb)
            return orbitals
        else:
            return NotImplemented
            
    def __neg__(self):
        orbitals = self.__class__()
        orbitals.extend(copy.deepcopy(self))
        for i,orb in enumerate(orbitals.copy()):
            coefs = [-coef for coef in orb.coefs]
            orbitals[i].coefs = coefs
        return orbitals

    def __mul__(self,other):
        if isinstance(other,(int,float)):
            orbitals = self.__class__()
            orbitals.extend(copy.deepcopy(self))
            for i,orb in enumerate(orbitals.copy()):
                coefs = [coef*other for coef in orb.coefs]
                orbitals[i].coefs = coefs
            return orbitals
        else:
            return NotImplemented

    def __rmul__(self,other):
        return(self.__mul__(other))

    def assign_symmetry(self,symmetry):
        """Assign symmetry to the orbital.

        Parameters
        ----------
        symmetry: Symmetry
            Symmetry of the orbital.

        Raises
        ------
        TypeError
            If symmetry is not Symmetry instance.
        """
        from moha.symmetry.symmetry import Symmetry
        if not isinstance(symmetry,Symmetry):
            raise TypeError("Parameters symmetry must be Symmetry instance.")
        if not symmetry.n in (0,1,2):
            raise ValueError("Number of electrons in orbital must be 0 1 or 2.")
        if not symmetry.ms2 in (0,1):
            raise ValueError("Spin polarization of orbital must be 0 or 1.")
        self.symmetry = symmetry
