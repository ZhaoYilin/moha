import numpy as np
import copy

__all__ = ['SlaterDeterminant']

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

