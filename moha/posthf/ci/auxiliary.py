from moha import *
import numpy as np
import copy
import itertools

class CI_Basis_Set(NElectronBasisSet):
    """
    """
    def __init__(self,hfwavefunction,truncation,size=0,basis_set=[]):
        super().__init__(size,basis_set)
        self.hfwavefunction = hfwavefunction
        self.dim = hfwavefunction.dim
        self.reference = hfwavefunction.configuration
        self.truncation = truncation

    def generate_basis_set(self,truncation):
        """Generate N electron basis set.

        Parameters:
            
        truncation
            type list
            a list of truncation degree
        """
        for i in truncation:
            list = self.generate_annihilation_creation(i)
            for j in list:
                confi = self.apply_annihilation_creation(j)
                slater = SlaterDeterminant(confi)
                self.add(slater)
        
        return self.basis_set

    def generate_annihilation_creation(self,degree):
        """
        """
        anni_list = []
        dim = self.dim
        
        occ_alpha_number = self.hfwavefunction.occ['alpha']
        occ_beta_number = self.hfwavefunction.occ['beta']
        
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
                            anni_list.append(dic) 

        return anni_list


    def apply_annihilation_creation(self,annihilation_creation):
        """
        """
        configuration = copy.deepcopy(self.reference)
        for i in annihilation_creation['annihilation_alpha']:
            configuration['alpha'][i] -= 1 

        for i in annihilation_creation['annihilation_beta']:
            configuration['beta'][i] -= 1

        for i in annihilation_creation['creation_alpha']:
            configuration['alpha'][i] += 1 

        for i in annihilation_creation['creation_beta']:
            configuration['beta'][i] += 1 
        return configuration



class CI_Hamiltonian(object):
    """
    """
    def __init__(self,ham,hfwavefunction,ci_basis_set):
        self.ham = ham
        self.hfwavefunction = hfwavefunction
        self.ci_basis_set = ci_basis_set
        self.h1e,self.h2e = self.molecular_spin_integral()

    def generate_matrix(self,ci_basis_set):
        """Generate CI Hamiltonian matrix.

        Parameters:

        ci_basis_set
            A CI_Basis_Set instance.
        """
        basis = ci_basis_set.basis_set
        num_det = ci_basis_set.size
        matrix = np.zeros((num_det,num_det))
        for i in range(num_det):
            for j in range(i+1):
                matrix[i,j] = self.calculate_matrix_element(basis[i],basis[j])
                matrix[j,i] = matrix[i,j]
        return matrix

    def calculate_matrix_element(self,det1,det2):
        """Calculate matrix element between two determinants.

        Parameters:

        det1
            A SlaterDeterminant instance.
        
        det2
            A SlaterDeterminant instance.
        """
        degree = det1.degree_of_excitation(det2)
        if degree == 0:
            value = self.calculate_matrix_element_identical(det1)
        elif degree == 1:
            value = self.calculate_matrix_element_diff1(det1,det2)
        elif degree == 2:
            value = self.calculate_matrix_element_diff2(det1,det2)
        else:
            value = 0.0
        return value

    def calculate_matrix_element_identical(self,det):
        """Calculate matrix element whose configurations are identical.

        Parameters:

        det
            A SlaterDeterminant instance.
        """
        O1 = 0.0
        O2 = 0.0
        
        occ_mixed_index = self.mixed_index_transformation(det.occ_index)
        for i in occ_mixed_index:
            O1 += self.h1e[i,i]

        for i in occ_mixed_index:
            for j in occ_mixed_index:
                O2 += self.h2e[i,j,i,j]

        value = O1+0.5*O2
        return value

    def calculate_matrix_element_diff1(self,det1,det2):
        """Calculate matrix element whose configuration differ by 1 spin orbitals.

        Parameters:

        det1
            A SlaterDeterminant instance.
        
        det2
            A SlaterDeterminant instance.
        """
        O1 = 0.0
        O2 = 0.0

        occ_mixed_index = self.mixed_index_transformation(det1.occ_index)
        sign = det1.phase(det2)
        anni = det1.annihilation_list(det2)
        crea = det1.creation_list(det2)

        a_list = self.mixed_index_transformation(anni)
        c_list = self.mixed_index_transformation(crea)
        for i in a_list:
            if i in occ_mixed_index:
                occ_mixed_index.remove(i)

        for j in c_list:
            if j in occ_mixed_index:
                occ_mixed_index.remove(j)

        O1 += self.h1e[a_list[0],c_list[0]]
        for k in occ_mixed_index:
            O2 += self.h2e[a_list[0],k,c_list[0],k]
        value = O1+O2
        return value*sign

    def calculate_matrix_element_diff2(self,det1,det2):
        """Calculate matrix element whose configuration differ by 2 spin orbitals.

        Parameters:

        det1
            A SlaterDeterminant instance.
        
        det2
            A SlaterDeterminant instance.
        """
        value = 0.0

        sign = det1.phase(det2)
        anni = det1.annihilation_list(det2)
        crea = det1.creation_list(det2)
        a_list = self.mixed_index_transformation(anni)
        c_list = self.mixed_index_transformation(crea)
        value += self.h2e[a_list[0],a_list[1],c_list[0],c_list[1]]
        return value*sign

    def mixed_index_transformation(self,dic):
        """Convert a configuration from dict form to list.

        Parameters:

        dic
            type dict
            A configuration in dict form
        """
        l = []
        for k in dic.keys():
            for i in dic[k]:
                if k=="alpha":
                    l.append(i*2)
                elif k=="beta":
                    l.append(i*2+1)
                else:
                    print("error")
        return l
        

    def molecular_spin_integral(self):
        """Compute the one and two electron integral in molecular spin orbital.

        output
            h1e h2e
        """
        ham = self.ham
        wf = self.hfwavefunction
        C = wf.coefficient
        ham.operators['kinetic'].basis_transformation(C)
        ham.operators['nuclear_attraction'].basis_transformation(C)
        h1e = ham.operators['kinetic'].spin
        h1e += ham.operators['nuclear_attraction'].spin
        ham.operators['electron_repulsion'].basis_transformation(C)
        h2e = ham.operators['electron_repulsion'].double_bar
        return h1e,h2e       
