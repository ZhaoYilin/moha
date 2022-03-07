from scipy.special import factorial2 as fact2

import numpy as np
import copy

class GaussianOrbital(object):
    """Gaussian type orbital class.
    
    Attributes
    ----------
    n_number : int
        Principal quantum number
    
    shell : list
        Angular momentum
    
    exp : scalar
        Primitive Gaussian exponent
    
    coef : scalar
        Primitive Gaussian coefficient
    
    norm : scalar
        Normalization factor
    
    origin : list
        Coordinate of the nuclei
    
    Methods
    -------

    """
    def __init__(self,type,atom_index,origin,n_number,shell=(),exps=[],coefs=[]):
        """Initialize the instance.

        Parameters
        ----------
        n_number : int
            Principal quantum number
    
        shell : list
            Angular momentum
    
        exp : scalar
            Primitive Gaussian exponent
    
        coef : scalar
            Primitive Gaussian coefficient
    
        norm : scalar
            Normalization factor
    
        origin : list
            Coordinate of the nuclei
        """
        self.type = type
        self.atom_index = atom_index
        self.origin = origin
        self.n_number = n_number
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
    def spatial(cls,atom_index,origin,n_number=0,shell=(),exps=[],coefs=[]):
        orbital = cls('spatial',atom_index,origin,n_number,shell,exps,coefs)
        return orbital

    @classmethod
    def spin(cls,atom_index,origin,spin,n_number=0,shell=(),exps=[],coefs=[]):
        orbital = cls('spin',atom_index,origin,n_number,shell,exps,coefs)
        orbital.spin = spin
        return orbital

