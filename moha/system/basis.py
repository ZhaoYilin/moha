import numpy as np
from .auxiliary import fact2
import copy

class OneElectronBasis(object):
    def __init__(self):
        pass

class NElectronBasis(object):
    def __init__(self):
        pass

class GaussianOrbital(OneElectronBasis):
    """
        n_number:
            type int
            principal quantum number
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
    def __init__(self,type,atom_index,origin,n_number,shell=(),exps=[],coefs=[]):
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
        type dict
        Binary representation of the determinants
        {'alpha':[1,1,0,0],'beta':[1,1,0,0]}
    """

    def __init__(self,configuration={}):
        """
        Constructor for the Determinant
        """
        self.configuration = configuration

    @property
    def configuration_alpha(self):
        """
        """
        return self.configuration['alpha']

    @property
    def configuration_beta(self):
        return self.configuration['beta']

    @property
    def occ_index(self):
        occ_alpha = []
        occ_beta = []
        for index,item in enumerate(self.configuration['alpha']):
            if item==1:
                occ_alpha.append(index)

        for index,item in enumerate(self.configuration['beta']):
            if item==1:
                occ_beta.append(index)

        return {"alpha":occ_alpha,"beta":occ_beta}

    def degree_of_excitation_alpha(self,det):
        """
        Compute the degree of excitation between D1 and D2 for alpha spin.
        """

        diff = list(np.array(det.configuration['alpha']) - np.array(self.configuration['alpha']))
        degree = int(sum(map(abs,diff))/2)
        return degree
    
    def degree_of_excitation_beta(self,det):
        """
        Compute the degree of excitation between D1 and D2 for beta spin.
        """

        diff = list(np.array(det.configuration['beta']) - np.array(self.configuration['beta']))
        degree = int(sum(map(abs,diff))/2)
        return degree

    def degree_of_excitation(self,det):
        """
        Compute the degree of excitation between D1 and D2.
        """
        alpha_degree = self.degree_of_excitation_alpha(det)
        beta_degree = self.degree_of_excitation_beta(det)
        degree = alpha_degree + beta_degree
        return degree

    def annihilation_list_alpha(self,det):
        """
        Identifying the annihilated alpha spin orbitals
        """
        diff = np.array(det.configuration['alpha']) - np.array(self.configuration['alpha'])
        if np.sum(np.abs(diff)) == 0: 
            index = []
        else:
            index = np.where(diff==-1)[0].tolist()
        return index


    def annihilation_list_beta(self,det):
        """
        Identifying the annihilated beta spin orbitals
        """
        diff = np.array(det.configuration['beta']) - np.array(self.configuration['beta'])
        if np.sum(np.abs(diff)) == 0: 
            index = []
        else:
            index = np.where(diff==-1)[0].tolist()
        return index

    def creation_list_alpha(self,det):
        """
        Identifying the created alpha spin orbitals
        """
        diff = np.array(det.configuration['alpha']) - np.array(self.configuration['alpha'])
        if np.sum(np.abs(diff)) == 0: 
            index = []
        else:
            index = np.where(diff==1)[0].tolist()
        return index

    def creation_list_beta(self,det):
        """
        Identifying the created beta spin orbitals
        """
        diff = np.array(det.configuration['beta']) - np.array(self.configuration['beta'])
        if np.sum(np.abs(diff)) == 0: 
            index = []
        else:
            index = np.where(diff==1)[0].tolist()
        return index

    def annihilation_list(self,det):
        """Identifying the annihilated spin orbitals
        det
            type object
            another Slater Determinant
        """
        alpha_list = self.annihilation_list_alpha(det)
        beta_list = self.annihilation_list_beta(det)
        index = {'alpha':alpha_list, 'beta':beta_list}
        return index
    
    def creation_list(self,det):
        """
        Identifying the created spin orbitals
        """
        alpha_list = self.creation_list_alpha(det)
        beta_list = self.creation_list_beta(det)
        index = {'alpha':alpha_list, 'beta':beta_list}
        return index

    def phase(self,det):
        """
        """
        tmp_configuration = copy.deepcopy(self.configuration)
        sign = 1
        anni_list = self.annihilation_list(det)
        crea_list = self.creation_list(det)
        for k in ("alpha","beta"):
            for i,j in zip(anni_list[k],crea_list[k]):
                if i <=j:
                    N_permutation = sum(tmp_configuration[k][i+1:j])
                else:    
                    N_permutation = sum(tmp_configuration[k][j+1:i])
                sign *= (-1)**(N_permutation)
                tmp_configuration[k][i]-=1
                tmp_configuration[k][j]+=1

        return sign

class NElectronBasisSet(object):
    """
    """
    def __init__(self,size=0,basis_set=[]):
        self.size = size
        self.basis_set = basis_set

    def add(self,NEbasis):
        self.size += 1
        self.basis_set.append(NEbasis)
