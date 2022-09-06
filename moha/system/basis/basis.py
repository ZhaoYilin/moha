import numpy as np
import copy

__all__ = ['Orbital','SlaterDeterminant']

class ContractedGaussian(object):
    """Predetermined linear combination of radial parts of GTOs
    Atomic orbtial represented by Contracted Gaussian functions.
    
    Attributes
    ----------
    origin : 3-list
        Coordinate of the nuclei.

    n : int
        Principal quantum number.
    
    shell : 3-tuple
        Angular momentum.
    
    coefs : list
        Primitive Gaussian coefficients.

    exps : list
        Primitive Gaussian exponents.
    
    Properties
    ----------
    norm: list
        Normalization factor.
    
    Methods
    -------
    __init__(self,type,origin,n,shell=(),coefs=[],exps=[])
        Initialize the instance.

    """
    def __init__(self,origin,n,shell=(),coefs=[],exps=[]):
        """Initialize the instance.

        Parameters
        ----------
        origin : list
            Coordinate of the nuclei.

        n : int
            Principal quantum number.
    
        shell : list
            Angular momentum.

        coefs : list
            Primitive Gaussian coefficients.
    
        exps : list
            Primitive Gaussian exponents.
        """
        self.origin = origin
        self.n = n
        self.shell = shell
        self.coefs = coefs
        self.exps  = exps

    def __neg__(self):
        coefs = [-coef for coef in self.coefs]
        return self.__class__(self.origin,self.n,self.shell,coefs,self.exps)
    
    def __add__(self,other):
        if self.parallel(other):
            coefs = [a+b for a,b in zip(self.coefs,other.coefs)]
            return self.__class__(self.origin,self.n,self.shell,coefs,self.exps)
        else:
            return NotImplemented 
        
    def __sub__(self,other):
        if self.parallel(other):
            coefs = [a-b for a,b in zip(self.coefs,other.coefs)]
            return self.__class__(self.origin,self.n,self.shell,coefs,self.exps)
        else:
            return NotImplemented 

    def __mul__(self,other):
        if isinstance(other,(int,float)):
            coefs = [coef*other for coef in self.coefs]
            return self.__class__(self.origin,self.n,self.shell,coefs,self.exps)
        elif isinstance(other,self.__class__):
            return NotImplemented
        else:
            return NotImplemented

    def __rmul__(self,other):
        return(self.__mul__(other))
    
    def __truediv__(self,other):
        if isinstance(other,(int,float)):
            coefs = [coef/other for coef in self.coefs]
            return self.__class__(self.origin,self.n,self.shell,coefs,self.exps)
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
        from scipy.special import factorial2 as fact2
        i,j,k = self.shell
        # self.norm is a list of length equal to number primitives
        norm = np.sqrt(np.power(2,2*(i+j+k)+1.5)*
                        np.power(self.exps,i+j+k+1.5)/
                        fact2(2*i-1)/fact2(2*j-1)/
                        fact2(2*k-1)/np.power(np.pi,1.5))
        return norm

    def parallel(self,other):
        if isinstance(other,self.__class__):
            self_attri = [self.origin,self.n,self.shell,self.exps]
            other_attri = [other.origin,other.n,other.shell,other.exps]
            if self_attri == other_attri:
                return True
            else:
                return False
        else:
            return False

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

class SlaterDeterminant(list):
    """Slater Determinant Basis.

    Attributes
    ----------
    self : list
        Binary representation of the determinant in spin orbitals.
        odd indices for alpha spin orbitals.
        even indices for beta spin orbitals.
        0 for unoccupied spin orbital.
        1 for occupied spin orbital.

    Properties
    ----------
    nspatials : int
        Number of spatial orbitals.
    
    nspins : int
        Number of spin orbitals.
    
    nalphas : int
        Number of alpha spin electrons.
    
    nbetas : int
        Number of beta spin electrons.
    
    nelectrons :
        Number of total electrons.
    
    alpha : list
        Configuraiton of alpha spin orbitals.
    
    beta : list
        Configuraiton of beta spin orbitals.
    
    occ_indices : list
        Occupation indices of the spin orbitals.
    
    spatial_format : list
        Representation of the Slater determinants in spatial orbitals.
    
    Methods
    -------
    __init__(self,configuration)
        Initialize the instance.
    
    assign_configuration(self,configuration)
        Assign configuration to the Slater determinant.
    
    degree_of_excitation(self,other)
        Compute the degree of excitation between D1 and D2.
    
    annihilation_list(self,other)
        Identifying the annihilated spin orbitals.
    
    creation_list(self,other)
        Identifying the created spin orbitals.
    
    phase(self,other)
        The phase is calculated as −1^{Nperm}.This number is equal
        to the number of occupied spin-orbitals between these two positions.
    """
    def __init__(self,configuration):
        """Initialize the instance.
        Parameters
        ----------
        configuration : list
            Representation of the Slater determinant.
        """
        self.assign_configuration(configuration)        

    def __repr__(self):
         return 'alpha: %s \n beta: %s' % (self.alpha, self.beta)
    
    @property
    def nspatials(self):
        """Number of spatial orbitals.
        """
        return int(len(self)/2)
    
    @property
    def nspins(self):
        """Number of spin orbitals.
        """
        return len(self)
    
    @property
    def nalphas(self):
        """Number of alpha spin electrons.
        """
        return sum(self.alpha)
    
    @property
    def nbetas(self):
        """Number of beta spin electrons.
        """
        return sum(self.beta)
    
    @property
    def nelectrons(self):
        """Number of total electrons.
        """
        return sum(self)
    
    @property
    def alpha(self):
        """Configuraiton of alpha spin orbitals.
        
        Returns
        -------
        alpha : list
            Configuraiton of alpha spin orbitals.
        """
        return self[0::2]

    @property
    def beta(self):
        """Configuraiton of beta spin orbitals.
        
        Returns
        -------
        beta : list
            Configuraiton of beta spin orbitals.
        """
        return self[1::2]
    
    @property
    def occ_indices(self):
        """Indices of occupied spin orbitals.
        
        Returns
        -------
        indices : list
            Indices of occupied spin orbitals.
        """
        indices = []
        for index,item in enumerate(self):
            if item==1:
                indices.append(index)
        return indices
    
    @property
    def vir_indices(self):
        """Indices of virtual spin orbitals.
        
        Returns
        -------
        vir_indices : list
            Indices of virtual spin orbitals.
        """
        indices = []
        for index,item in enumerate(self):
            if item==0:
                indices.append(index)
        return indices
    
    @property
    def spin_format(self):
        """Representation of the Slater determinants in spin orbitals.
        
        Returns
        -------
        spin_format : list
            Representation of the Slater determinants in spin orbitals.
        """
        return list(self)
    
    @property
    def spatial_format(self):
        """Representation of the Slater determinants in spatial orbitals.
            [0,0,0,2,1,3,3,3]
            0 for unoccupied orbital
            1 for occupied alpha spin orbital
            2 for occupied beta spin orbital
            3 for doubly occupied orbital
        
        Returns
        -------
        spatial_format : list
            Representation of the Slater determinants in spatial orbitals.
        """
        return list(np.array(self.alpha)+ 2*np.array(self.beta))
    
    @classmethod
    def ground(cls, nocc, nspins):
        """Create a ground state Slater determinant.
        
        Parameters
        ----------
        nocc : int
            Number of occupied spin orbitals.

        norbs : int
            Total number of spin orbitals.
        
        Returns
        -------
        ground_sd : gmpy2.mpz
            Integer that describes the occupation of a Slater determinant as a
            bitstring.
        
        Raises
        ------
        ValueError
            If `nocc > nspins`.
        
        Notes
        -----
        Assumes that the spin-orbitals are ordered by energy from lowest to
        greatest. The occupation is assumed to have the alpha block frist, 
        then the beta block.
        """
        from math import floor, ceil
        nspatials = int(nspins/2)
        alpha = [1]*ceil(nocc/2) + [0]*(nspatials-ceil(nocc/2))
        beta = [1]*floor(nocc/2) + [0]*(nspatials-floor(nocc/2))
        configuration = []
        for i,j in zip(alpha, beta):
            configuration.append(i)
            configuration.append(j)

        return cls(configuration)
    
    def assign_configuration(self,configuration):
        """Assign configuration to the Slater determinant.
        
        Parameters
        ----------
        configuration : list
            Binary representation of the Slater determinant.
        
        Raises
        ------
        TypeError
            If configuration is not list.
        ValueError
            If configuration is not list.
        """
        if not isinstance(configuration,list):
            raise TypeError("Configuration must be list.")
        if not len(configuration)%2==0:
            raise TypeError("Length of the configuration in spin orbitals must be even.")
        list.__init__(self,configuration)

    def degree_of_excitation(self,other):
        """Compute the degree of excitation between D1 and D2 for alpha spin.

        Parameters
        ----------
        other : SlaterDeterminant
            The second Slater determinant instance.
        
        Raises
        ------
        TypeError
            If det is not a SlaterDeterminant instance.

        Returns
        -------
        degree : int
            Degree of alpha spin excitation between two Slater determinant.
        """
        if not isinstance(other,SlaterDeterminant):
            raise TypeError("Parameter other must be a SlaterDeterminant instance.")
        diff = list(np.array(other) - np.array(self))
        degree = int(sum(map(abs,diff))/2)
        return degree

    def annihilation_list(self,other):
        """Identifying the annihilated alpha spin orbitals.

        Parameters
        ----------
        other : SlaterDeterminant
            The second Slater determinant instance.
        
        Raises
        ------
        TypeError
            If det is not a SlaterDeterminant instance.

        Returns
        -------
        index : dict
            Index of annihilation alpha spin.
        """
        if not isinstance(other,SlaterDeterminant):
            raise TypeError("Parameter other must be a SlaterDeterminant instance.")
        diff = np.array(other) - np.array(self)
        if np.sum(np.abs(diff)) == 0: 
            indices = []
        else:
            indices = np.where(diff==-1)[0].tolist()
        return indices

    def creation_list(self,other):
        """Identifying the created alpha spin orbitals.
        
        Parameters
        ----------
        other : SlaterDeterminant
            The second Slater determinant instance.
        
        Raises
        ------
        TypeError
            If det is not a SlaterDeterminant instance.

        Returns
        -------
        index : dict
            Index of creation alpha spin.
        """
        if not isinstance(other,SlaterDeterminant):
            raise TypeError("Parameter other must be a SlaterDeterminant instance.")
        diff = np.array(other) - np.array(self)
        if np.sum(np.abs(diff)) == 0: 
            indices = []
        else:
            indices = np.where(diff==1)[0].tolist()
        return indices

    def phase(self,other):
        """The phase is calculated as −1^{Nperm}.This number is equal
        to the number of occupied spin-orbitals between these two positions.
        
        Parameters
        ----------
        other : SlaterDeterminant
            The second Slater determinant instance.
        
        Raises
        ------
        TypeError
            If det is not a SlaterDeterminant instance.

        Returns
        -------
        sign : int
            Phase for the evaluation between two Slater determinant.
        """
        if not isinstance(other,SlaterDeterminant):
            raise TypeError("Parameter other must be a SlaterDeterminant instance.")
        tmp_configuration = copy.deepcopy(self)
        sign = 1
        anni_list = self.annihilation_list(other)
        crea_list = self.creation_list(other)
        for i,j in zip(anni_list,crea_list):
            if i <=j:
                N_permutation = sum(tmp_configuration[i+1:j])
            else:    
                N_permutation = sum(tmp_configuration[j+1:i])
            sign *= (-1)**(N_permutation)
            tmp_configuration[i]-=1
            tmp_configuration[j]+=1
        return sign

