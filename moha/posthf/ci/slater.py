import numpy as np
import copy

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
        return list(map(add, self.alpha, 2*self.beta))
    
    @classmethod
    def ground(cls, nocc, nspins):
        """Create a ground state Slater determinant.
        If the number of electrons is odd, then the last electron is put into an alpha orbital.
        
        Parameters
        ----------
        nocc : int
            Number of occupied spin orbitals.
        norbs : int
            Total number of spin orbitals.
        
        Returns
        -------
        ground_sd : gmpy2.mpz
            Integer that describes the occupation of a Slater determinant as a bitstring.
        
        Raises
        ------
        ValueError
            If `nocc > nspins`.
        
        Notes
        -----
        Assumes that the spin-orbitals are ordered by energy from lowest to greatest. The occupation is
        assumed to have the alpha block frist, then the beta block.
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
        det : SlaterDeterminant
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
        det : SlaterDeterminant
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
        det : SlaterDeterminant
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

class SlaterCondon(object):
    """The rules which allow us to express matrix elements Hij in terms of one and two
    electron integrals.
    
    Attributes
    ----------
    h1e : np.ndarray
        One electron integral.

    g2e : np.ndarray
        Two electron integral.
    
    Methods
    -------
    __init__(self, h1e, g2e)
        Initialize the instance.

    calculate_matrix_element(self,det1,det2)
        Calculate matrix element between two determinants.
    
    identical(self,det)
        Calculate matrix element whose configurations are identical.
    
    differ_by1(self,det1,det2):
        Calculate matrix element whose configuration differ by 1 spin orbitals.
    
    difer_byf2(self,det1,det2):
        Calculate matrix element whose configuration differ by 2 spin orbitals.
    """
    def __init__(self, h1e, g2e):
        """
        Parameters
        ----------
        h1e : np.ndarray
            One electron integral.

        g2e : np.ndarray
            Two electron integral.
        """
        self.h1e = h1e
        self.g2e = g2e
    
    def calculate_matrix_element(self,sd1,sd2):
        """Calculate matrix element between two determinants.

        Parameters
        ----------
        sd1 : SlaterDeterminant
            The first SlaterDeterminant instance.
        
        sd2 : SlaterDeterminant
            The second SlaterDeterminant instance.
        
        Raises
        ------
        TypeError
            If sd1 is not SlaterDeterminant instance.
            If sd2 is not SlaterDeterminant instance.

        Returns
        -------
        value : float
            The integral value of CI Hamiltonian.
        """
        if not isinstance(sd1,SlaterDeterminant):
            raise TypeError("Parameters det1 must be SlaterDeterminant instance.")
        if not isinstance(sd2,SlaterDeterminant):
            raise TypeError("Parameters det2 must be SlaterDeterminant instance.")
        degree = sd1.degree_of_excitation(sd2)
        if degree == 0:
            value = self.identical(sd1)
        elif degree == 1:
            value = self.differ_by1(sd1,sd2)
        elif degree == 2:
            value = self.differ_by2(sd1,sd2)
        else:
            value = 0.0
        return value
    
    def identical(self,sd):
        """Calculate matrix element whose configurations are identical.

        Parameters
        ----------
        sd : SlaterDeterminant
            A SlaterDeterminant instance.
        
        Raises
        ------
        TypeError
            If sd is not SlaterDeterminant instance.
        
        Returns
        -------
        value : float
            The integral value of CI Hamiltonian.
        """
        if not isinstance(sd, SlaterDeterminant):
            raise TypeError("Parameters sd must be SlaterDeterminant instance.")
        O1 = 0.0
        O2 = 0.0
        for i in sd.occ_indices:
            O1 += self.h1e[i,i]
        for i in sd.occ_indices:
            for j in sd.occ_indices:
                O2 += self.g2e[i,j,i,j]
        value = O1+0.5*O2
        return value
    
    def differ_by1(self,sd1,sd2):
        """Calculate matrix element whose configuration differ by 1 spin orbitals.

        Parameters
        ----------
        sd1
            A SlaterDeterminant instance.
        
        sd2
            A SlaterDeterminant instance.
        
        Raises
        ------
        TypeError
            If sd1 is not SlaterDeterminant instance.
            If sd2 is not SlaterDeterminant instance.
        
        Returns
        -------
        value : float
            The integral value of CI Hamiltonian.
        """
        if not isinstance(sd1,SlaterDeterminant):
            raise TypeError("Parameters sd1 must be SlaterDeterminant instance.")
        if not isinstance(sd2,SlaterDeterminant):
            raise TypeError("Parameters sd2 must be SlaterDeterminant instance.")
        O1 = 0.0
        O2 = 0.0

        sign = sd1.phase(sd2)
        occ_indices = copy.deepcopy(sd1.occ_indices)
        a_list = sd1.annihilation_list(sd2)
        c_list = sd1.creation_list(sd2)

        for i in a_list:
            if i in occ_indices:
                occ_indices.remove(i)

        for j in c_list:
            if j in occ_indices:
                occ_indices.remove(j)
        O1 += self.h1e[a_list[0],c_list[0]]
        for k in occ_indices:
            O2 += self.g2e[a_list[0],k,c_list[0],k]
        value = O1+O2
        value *=sign
        return value

    def differ_by2(self,sd1,sd2):
        """Calculate matrix element whose configuration differ by 2 spin orbitals.

        Parameters
        ----------
        sd1
            First SlaterDeterminant instance.
        
        sd2
            Second SlaterDeterminant instance.
        
        Raises
        ------
        TypeError
            If sd1 is not SlaterDeterminant instance.
            If sd2 is not SlaterDeterminant instance.

        Returns
        -------
        value : float
            The integral value of CI Hamiltonian.
        """
        if not isinstance(sd1,SlaterDeterminant):
            raise TypeError("Parameters sd1 must be SlaterDeterminant instance.")
        if not isinstance(sd2,SlaterDeterminant):
            raise TypeError("Parameters sd2 must be SlaterDeterminant instance.")
        value = 0.0

        sign = sd1.phase(sd2)
        a_list = sd1.annihilation_list(sd2)
        c_list = sd1.creation_list(sd2)
        value += self.g2e[a_list[0],a_list[1],c_list[0],c_list[1]]
        value *=sign
        return value


