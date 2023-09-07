import itertools
import math
import numpy as np
from scipy.special import factorial2 as fact2

__all__ = ['PrimitiveGaussian','ContractedGaussian']

class PrimitiveGaussian(object):
    """Primitive Gaussian functions.
    
    Attributes
    ----------
    coefficient : float
        Primitive Gaussian coefficient.

    origin : List[float,float,float]
        Coordinate of the nuclei.

    shell : List[int,int,int]
        Angular momentum.

    exponent : float
        Primitive Gaussian exponent.
    
    Properties
    ----------
    norm: float
        Normalization factor.
    
    Methods
    -------
    __init__(self,type,origin,n,shell,coef,exp)
        Initialize the instance.

    __call__(self,x,y,z)
        Call the instance as function.
    """
    def __init__(self,coefficient,origin,shell,exponent):
        """Initialize the instance.

        Parameters
        ----------
        coefficient : float
            Primitive Gaussian coefficient.

        origin : List[float,float,float]
            Coordinate of the nuclei.

        shell : List[int,int,int]
            Angular momentum.

        exponent : float
            Primitive Gaussian exponent.
        """
        self.coefficient = coefficient
        self.origin = origin
        self.shell = shell
        self.exponent  = exponent
        self.norm_memo = None

    def __call__(self,x,y,z):
        """Call the instance as function.
        """
        X = x-self.origin[0]
        Y = y-self.origin[1]
        Z = z-self.origin[2]
        rr = X**2+Y**2+Z**2
        result = np.power(X,self.shell[0])*\
            np.power(Y,self.shell[1])*\
            np.power(Z,self.shell[2])*\
            np.exp(-self.exponent*rr)
        return result

    @property
    def norm(self):
        """Normalization factors. 
        
        Return
        ------
        norm : list
            Normalization factors
        """
        if self.norm_memo is None:
            i,j,k = self.shell
            term1 = (2*self.exponent/math.pi)**(3/2)
            term2 = (4*self.exponent)**(i+j+k)
            term3 = fact2(2*i-1)*fact2(2*j-1)*fact2(2*k-1)
            self.norm_memo = math.sqrt(term1*term2/term3)
        return self.norm_memo

    def derivative(self):
        """ Differentiation of Cartesian Gaussians.
        """
        dx = []
        dy = []
        dz = []

        c = self.coefficient
        origin = self.origin
        l = self.shell[0]
        m = self.shell[1]
        n = self.shell[2]
        a = self.exponent
        
        # Dx
        primitive_basis = self.__class__(2*a*c,origin,(l+1,m,n),a)
        dx.append(primitive_basis)
        if l >=1:
            primitive_basis = self.__class__(-l*c,origin,(l-1,m,n),a)
            dx.append(primitive_basis)

        # Dy
        primitive_basis = self.__class__(2*a*c,origin,(l,m+1,n),a)
        dy.append(primitive_basis)
        if m >=1:
            primitive_basis = self.__class__(-m*c,origin,(l,m-1,n),a)
            dy.append(primitive_basis)
        
        # Dz
        primitive_basis = self.__class__(2*a*c,origin,(l,m,n+1),a)
        dz.append(primitive_basis)
        if n >=1:
            primitive_basis = self.__class__(-n*c,origin,(l,m,n-1),a)
            dz.append(primitive_basis)

        return dx, dy, dz

class ContractedGaussian(object):
    """Predetermined linear combination of radial parts of GTOs
    Atomic orbtial represented by Contracted Gaussian functions.
    
    Attributes
    ----------
    coefficients : List of float
        Primitive Gaussian coefficient.

    origin : List[float,float,float]
        Coordinate of the nuclei.

    shell : List[int,int,int]
        Angular momentum.

    exponents : List of float
        Primitive Gaussian exponent.
    
    Properties
    ----------
    norm: list
        Normalization factor.
    
    Methods
    -------
    __init__(self,type,origin,n,shell=(),coefs=[],exps=[])
        Initialize the instance.

    """
    def __init__(self,coefficients=[],origin=[],shell=[],exponents=[]):
        """Initialize the instance.

        Parameters
        ----------
        coefficients : List of float
            Primitive Gaussian coefficient.

        origin : List[float,float,float]
            Coordinate of the nuclei.

        shell : List[int,int,int]
            Angular momentum.

        exponents : List of float
            Primitive Gaussian exponent.
        """
        self.coefficients = coefficients
        self.origin = origin
        self.shell = shell
        self.exponents  = exponents
        self.norm_memo = None

    def __call__(self,x,y,z):
        X = x-self.origin[0]
        Y = y-self.origin[1]
        Z = z-self.origin[2]
        rr = X**2+Y**2+Z**2

        cg = np.zeros(x.shape)
        for coef, exp in zip(self.coefficients,self.exponents):
            cg += np.power(X,self.shell[0])*\
                np.power(Y,self.shell[1])*\
                np.power(Z,self.shell[2])*\
                np.exp(-exp*rr)
        return cg

    def __neg__(self):
        coefficients = [-coefficient for coefficient in self.coefficients]
        origin = self.origin
        shell = self.shell
        exponents = self.exponents
        new_CGTO = self.__class__(coefficients,origin,shell,exponents)
        return new_CGTO

    def __add__(self,other):
        if self.parallel(other):
            coefficients = [a+b for a,b in zip(self.coefficients,other.coefficients)]
            origin = self.origin
            shell = self.shell
            exponents = self.exponents
            new_CGTO = self.__class__(coefficients,origin,shell,exponents)
            return new_CGTO
        else:
            return NotImplemented

    def __sub__(self,other):
        if self.parallel(other):
            coefficients = [a-b for a,b in zip(self.coefficients,other.coefficients)]
            origin = self.origin
            shell = self.shell
            exponents = self.exponents
            new_CGTO = self.__class__(coefficients,origin,shell,exponents)
            return new_CGTO
        else:
            return NotImplemented

    def __mul__(self,other):
        if isinstance(other,(int,float)):
            coefficients = [other*coefficient for coefficient in self.coefficients]
            origin = self.origin
            shell = self.shell
            exponents = self.exponents
            new_CGTO = self.__class__(coefficients,origin,shell,exponents)
            return new_CGTO
        elif isinstance(other,self.__class__):
            return NotImplemented
        else:
            return NotImplemented

    def __rmul__(self,other):
        return(self.__mul__(other))

    def __truediv__(self,other):
        if isinstance(other,(int,float)):
            coefficients = [coefficient/other for coefficient in self.coefficients]
            origin = self.origin
            shell = self.shell
            exponents = self.exponents
            new_CGTO = self.__class__(coefficients,origin,shell,exponents)
            return new_CGTO
        elif isinstance(other,self.__class__):
            return NotImplemented
        else:
            return NotImplemented

    @property
    def norm(self):
        """Normalization factors. 
        
        Return
        ------
        norm : list
            Normalization factors
        """
        if self.norm_memo is None:
            l, m, n = self.shell

            result = 0.0
            for pga, pgb in itertools.product(self.expansion, repeat=2):
                a = pga.exponent
                b = pgb.exponent
                ca = pga.coefficient
                cb = pgb.coefficient
                na = pga.norm
                nb = pgb.norm

                term1 = fact2(2 * l - 1) * fact2(2 * m - 1) * fact2(2 * n - 1)
                term2 = (math.pi / (a + b)) ** (3 / 2)
                term3 = (2 * (a + b)) ** (l + m + n)
                result += (ca * cb * na * nb * term1 * term2) / term3
            self.norm_memo = 1 / math.sqrt(result)

        return self.norm_memo

    @property
    def expansion(self):
        """Normalization factors. 
        
        Return
        ------
        primitives : list
            Normalization factors
        """
        n = len(self.coefficients)
        primitives = []
        for i in range(n):
            coefficient = self.coefficients[i]
            exponent = self.exponents[i]
            pg = PrimitiveGaussian(coefficient,self.origin,self.shell,exponent)
            primitives.append(pg)
        return primitives

    def parallel(self,other):
        if isinstance(other,self.__class__):
            attributes_self = [self.origin,self.shell,self.exponents]
            attributes_other = [other.origin,other.shell,other.exponents]
            if attributes_self == attributes_other:
                return True
            else:
                return False
        else:
            return False

