from moha.system.basis.basis import Orbital,SlaterDeterminant
from moha.io.log import log

import copy
import itertools

__all__ = ['BasisSet','FockSpace']

class BasisSet(list):
    """One electron basis set.
    
    Classmethods
    ------------
    build(cls,molecule,basisfile)
        Build the instance.

    Methods
    -------
    __new__(cls)
        Generate new instance.
    
    __init__(self)
        Initialize the instance.
    
    append(self,basis)
        Apend basis to basis set.
    """
    def __new__(cls):
        """Generate new instance.
        """
        obj = list.__new__(cls)
        return obj

    @classmethod
    def build(cls,molecule,basisfile):
        """Build the basis set.

        Parameters
        ----------
        molecule: Molecule
            Moleucle instance.

        basisfile : str
            Name of the basis set file.

        Return
        ------
        basis_set: BasisSet
            Basis set instance.
        """
        log.hline(char='=')
        log.blank()
        log('Basis Set Section'.format())
        log.hline()

        basis_set = BasisSet()

        log('{0:<20s}{1:<20s}{2:<20s}'.format('Element','Number','Type'))
        log.hline()
        if not isinstance(basisfile, str):
            raise TypeError("basisfile must be str")
        if basisfile.endswith('.nwchem'):
            from moha.system.basis import nwchem
            for atom in molecule:
                old_len = len(basis_set)
                nwchem.load(atom,basis_set,basisfile)
                log('{0:<20s}{1:<20d}{2:<20s}'.format(atom.symbol,len(basis_set)-old_len,
                basisfile.strip('.nwchem')))
        else:
            raise ValueError("unknown file format")

        log.blank()
        return basis_set

    def append(self,basis):
        """Append basis to basis set.

        Parameters
        ----------
        basis : Orbital
            Orbital instance.
        
        Raises
        ------
        TypeError
            If basis is not an Orbital instance.
        """
        if not isinstance(basis, Orbital):
            raise TypeError("Basis must be Orbital instance")
        super().append(basis)

class FockSpace(list):
    """Fock space.

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
    
    _generate_annihilation_creation(self,degree)
        Generate a list of annihilation and creation.
    
    _apply_annihilation_creation(self,ac_list)
        Apply annihilation and creation to the self Slater determinant.
    """
    def __new__(cls):
        """Generate new basis set.
        """
        obj = super().__new__(cls)
        return obj

    @classmethod
    def subspace(cls,m,n):
        """The Fock subspace F(M,N).

        Parameters
        ----------
        m: int
            Number of spin orbitals.
        
        n: int
            Numer of electrons.

        Return
        ------
        fock_space: FockSpace
            The Fock subsapce F(M,N).
        """
        fock_space = cls()
        return fock_space

    @classmethod
    def ci_space(cls, reference, truncation):
        """Generate CI basis set.

        Parameters
        ----------
        truncation : int
            Truncation number.

        Return
        ------
        fock_space: FockSpace
            The Fock subsapce.

        Raises
        ------
        TypeError 
            If truncation is not a list.
        """
        fock_space = cls()

        if not isinstance(truncation, int):
            raise TypeError("Parameter truncation must be int.")
        for nexc in range(truncation+1):
            fock_space.add_excitation(reference,nexc)
            print(len(fock_space))
        return fock_space

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

        ac_list = self._generate_annihilation_creation(reference, nexc)
        for j in ac_list:
            print(j)
            confi = self._apply_annihilation_creation(reference, j)
            slater = SlaterDeterminant(confi)
            self.append(slater)

    def _generate_annihilation_creation(self, reference, degree):
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

    def _apply_annihilation_creation(self, reference, ac):
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