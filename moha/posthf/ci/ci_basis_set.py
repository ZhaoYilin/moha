from moha.posthf.ci.slater import SlaterDeterminant

import numpy as np
import copy
import itertools

class CIBasisSet(list):
    """Configuration interaction basis set.

    Attributes
    ----------
    reference
        Reference Slater Determinant.
        
    truncation : int
        Maxium excitation level.
    
    Methods
    -------
     __init__(self)
        Initialize the instance.
    
    append(self,basis)
        Append one basis to the basis set.    
    
    assign_truncation(self,truncation)
        Assign truncation.
    
    generate_basis_set(self,truncation)
        Generate N electron basis set.
    
    generate_annihilation_creation(self,degree)
        Generate a list of annihilation and creation.
    
    apply_annihilation_creation(self,ac_list)
        Apply annihilation and creation to the self Slater determinant.
    """
    def __new__(cls):
        """Generate new basis set.
        """
        obj = super().__new__(cls)
        
        return obj
    
    def __init__(self):
        """Initialize the CI basis set.

        Parameters
        ----------
        reference : SlaterDeterminant
            Reference configuration.
        """
        pass
    
    @property
    def size(self):
        """Size of the CI basis set.
        """
        return len(self)
    
    @classmethod
    def build(cls, reference, truncation):
        """Generate CI basis set.

        Parameters
        ----------
        truncation : int
            Truncation number.

        Raises
        ------
        TypeError 
            If truncation is not a list.
        """
        basis_set = cls()

        if not isinstance(truncation, int):
            raise TypeError("Parameter truncation must be int.")
        for nexc in range(truncation+1):
            basis_set.add_excitation(reference,nexc)
        
        return basis_set

    def append(self,basis):
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
        super().append(basis)  
    
    def add_excitation(self,reference,nexc):
        """Add excitation state.

        Parameters
        ----------
        reference : SlaterDeterminant
            Reference configuration.
        
        nexc : int
            Number of excitation electrons.
        
        Raises
        ------
        TypeError
            If reference is not a SlaterDeterminant instance.
        
        TypeError 
            If nexc is not int.
        """
        if not isinstance(reference,SlaterDeterminant):
            raise TypeError("Reference must be SlaterDeterminant instance.")
        if not isinstance(nexc,int):
            raise TypeError("Number of excitation electron must be int.")
        
        ac_list = self.generate_annihilation_creation(reference, nexc)
        for j in ac_list:
            confi = self.apply_annihilation_creation(reference, j)
            slater = SlaterDeterminant(confi)
            self.append(slater)
    
    def generate_annihilation_creation(self, reference, degree):
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
        nspatials = reference.nspatials
        occ_indices = reference.occ_indices
        vir_indices = reference.vir_indices

        for alpha_degree in range(degree+1):
            beta_degree = degree - alpha_degree
            annihilation_index_alpha = list(itertools.combinations(occ_indices[:nspatials],alpha_degree))
            annihilation_index_beta = list(itertools.combinations(occ_indices[nspatials:],beta_degree))
            creation_index_alpha = list(itertools.combinations(vir_indices[:nspatials],alpha_degree))
            creation_index_beta = list(itertools.combinations(vir_indices[nspatials:],beta_degree))
            for i in annihilation_index_alpha:
                for j in annihilation_index_beta:
                    for k in creation_index_alpha:
                        for l in creation_index_beta:
                            ac_list.append([list(i+j),list(k+l)])
        return ac_list


    def apply_annihilation_creation(self, reference, ac):
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
        if not isinstance(ac,list):
            raise TypeError("Parameter ac must be list.")
        configuration = copy.deepcopy(reference)
        
        for i in ac[0]:
            configuration[i] -= 1 
        for i in ac[1]:
            configuration[i] += 1 
        return configuration

