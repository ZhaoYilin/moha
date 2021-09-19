from moha.system.basis.slater_determinant import SlaterDeterminant
from moha.system.wavefunction.hf_wavefunction import HFWaveFunction
from moha.system.basis_set.base import BaseBasisSet
import numpy as np
import copy
import itertools

class CIBasisSet(BaseBasisSet):
    """Configuration interaction basis set.

    Attributes
    ----------
    reference
        Reference Hartree Fock wavefunction.
        
    truncation : list
        A list of truncation level.
    
    size : int
        Size of the basis set.

    bases : list
        A list of basis.

    Methods
    -------
     __init__(self,reference,truncation,size=0,bases=[])
        Initialize the instance.
    
    assign_size(self,size)
        Assign size to the basis set.

    assign_bases(self,bases)
        Assign basis to the basis set.
    
    append_basis(self,basis)
        Append one basis to the basis set.    
    
    assign_reference(self,reference)
        Assign reference Hartree Fock wavefunction.
    
    assign_truncation(self,truncation)
        Assign truncation.
    
    generate_basis_set(self,truncation)
        Generate N electron basis set.
    
    generate_annihilation_creation(self,degree)
        Generate a list of annihilation and creation.
    
    apply_annihilation_creation(self,ac_list)
        Apply annihilation and creation to the self Slater determinant.
    """
    def __init__(self,reference,truncation,size=0,bases=[]):
        """Initialize the configuration interaction basis set.

        Parameters
        ----------
        reference
            Reference Hartree Fock wavefunction.
        
        truncation : list
            A list of truncation level.
        
        size : int
            Size of the basis set.

        bases : list
            A list of the basis.
        """
        super().__init__(size,[])
        self.assign_reference(reference)
        self.assign_truncation(truncation)
        self.generate_basis_set(truncation)

    def append_basis(self,basis):
        """Add one basis to the basis set.

        Parameters
        ----------
        basis
            Basis instance.

        Raises
        ------
        TypeError
            If orb is not a SlaterDeterminant instance.
        """
        if not isinstance(basis, SlaterDeterminant):
            raise TypeError("Basis must be a SlaterDeterminant instance")
        self.size += 1
        self.bases.append(basis)  

    def assign_reference(self,reference):
        """Assign reference Hartree Fock wavefunction.

        Parameters
        ----------
        reference
            Reference Hartree Fock wavefunction.
        
        Raises
        ------
        TypeError
            If reference is not a Hartree Fock wavefunction instance.
        """
        if not isinstance(reference,HFWaveFunction):
            raise TypeError("Reference must be a HFWaveFunction instance.")        
        self.reference = reference
    
    def assign_truncation(self,truncation):
        """Assign truncation.

        Parameters
        ----------
        truncation : list
            A list of truncation level.
        
        Raises
        ------
        TypeError 
            If truncation is not a list.
        """
        if not isinstance(truncation,list):
            raise TypeError("Parameter truncation must be a list.")
        self.truncation = truncation

    def generate_basis_set(self,truncation):
        """Generate N electron basis set.

        Parameters
        ----------
        truncation : list
            A list of truncation degree.

        Raises
        ------
        TypeError 
            If truncation is not a list.
        """
        if not isinstance(truncation,list):
            raise TypeError("Parameter truncation must be a list.")
        for i in truncation:
            ac_list = self.generate_annihilation_creation(i)
            for j in ac_list:
                confi = self.apply_annihilation_creation(j)
                slater = SlaterDeterminant(confi)
                self.append_basis(slater)

    def generate_annihilation_creation(self,degree):
        """Generate a list of annihilation and creation.

        Parameters
        ----------
        degree : int
            Truncation level.
        
        Raises
        ------
        TypeError
            If parameter degree is not an int.

        ValueError
            If parameter degree is not a none negative number.
            
        Returns
        -------
        ac_list : list
            A list of annihilation and creation.
        """
        if not isinstance(degree,int):
            raise TypeError("Parameter degree must be an int.")
        if degree<0:
            raise ValueError("Parameter degree must be a none negative number.")
        ac_list = []
        dim = self.reference.nspatial
        
        occ_alpha_number = self.reference.occ['alpha']
        occ_beta_number = self.reference.occ['beta']
        
        occ_alpha_index = list(range(occ_alpha_number))
        occ_beta_index = list(range(occ_beta_number))
        
        uocc_alpha_index = list(range(occ_alpha_number,dim))
        uocc_beta_index = list(range(occ_beta_number,dim))
        
        for alpha_degree in range(degree+1):
            annihilation_index_alpha = list(itertools.combinations(occ_alpha_index,alpha_degree))
            creation_index_alpha = list(itertools.combinations(uocc_alpha_index,alpha_degree))
            beta_degree = degree - alpha_degree
            annihilation_index_beta = list(itertools.combinations(occ_beta_index,beta_degree))
            creation_index_beta = list(itertools.combinations(uocc_beta_index,beta_degree))
            for i in annihilation_index_alpha:
                for j in annihilation_index_beta:
                    for k in creation_index_alpha:
                        for l in creation_index_beta:
                            dic = {"annihilation_alpha":list(i),"annihilation_beta":list(j)
                            ,"creation_alpha":list(k),"creation_beta":list(l)}
                            ac_list.append(dic) 
        return ac_list


    def apply_annihilation_creation(self,ac):
        """Apply annihilation and creation to the self Slater determinant.

        Parameters
        ----------
        ac : dict
            A dictionary of annihilation and creation.
        
        Raises
        ------
        TypeError
            If ac is not a dictionary.

        Returns
        -------
        configuraiton : dict
            New configuration for the self Slater determinant.
        """
        if not isinstance(ac,dict):
            raise TypeError("Parameter ac must be a dictionary.")
        configuration = copy.deepcopy(self.reference.configuration)
        for i in ac['annihilation_alpha']:
            configuration['alpha'][i] -= 1 
        for i in ac['annihilation_beta']:
            configuration['beta'][i] -= 1
        for i in ac['creation_alpha']:
            configuration['alpha'][i] += 1 
        for i in ac['creation_beta']:
            configuration['beta'][i] += 1 
        return configuration

