import numpy as np
from .auxiliary import fact2
class OneElectronBasis(object):
    def __init__(self):
        pass

class NElectronBasis(object):
    def __init__(self):
        pass

class GaussianOrbital(OneElectronBasis):
    """
        shell: 
            type list
            angular momentum
        exp:   
            type scalar
            primitive Gaussian exponent
        coef:  
            type scalar
            primitive Gaussian coefficient
        norm:   
            type scalar
            normalization factor
        origin:
            type list
            coordinate of the nuclei
    """
    def __init__(self,type,atom_index,origin,shell=(),exps=[],coefs=[]):
        self.type = type
        self.atom_index = atom_index
        self.origin = origin
        self.shell = shell
        self.exps  = exps
        self.coefs = coefs
        self.norm = self.normalize()

    def normalize(self):
        """ method to calculate the normalization factors 
        """
        l,m,n = self.shell
        # self.norm is a list of length equal to number primitives
        norm = np.sqrt(np.power(2,2*(l+m+n)+1.5)*
                        np.power(self.exps,l+m+n+1.5)/
                        fact2(2*l-1)/fact2(2*m-1)/
                        fact2(2*n-1)/np.power(np.pi,1.5))
        return norm

    @classmethod
    def spatial(cls,atom_index,origin,shell=(),exps=[],coefs=[]):
        orbital = cls('spatial',atom_index,origin,shell,exps,coefs)
        return orbital

    @classmethod
    def spin(cls,atom_index,origin,spin,shell=(),exps=[],coefs=[]):
        orbital = cls('spin',atom_index,origin,shell,exps,coefs)
        orbital.spin = spin
        return orbital

class BasisSet(object):
    def __init__(self,size=0,basis=[]):
        self.size = size
        self.basis = basis

    def add(self,orb):
        self.size += self.size
        self.basis.append(orb)

class SlaterDeterminant(NElectronBasis):
    """
    configuration:
        {'alpha':[1,1,0,0],'beta':[1,1,0,0]}
    """
    def __init__(self,configuration):
        self.configuration = configuration

class NElectronBasisSet(object):
    """
    """
    def __init__(self,size=0,basis=[]):
        self.size = size
        self.basis = basis

    def add(self,orb):
        self.size += 1
        self.basis.append(orb)
