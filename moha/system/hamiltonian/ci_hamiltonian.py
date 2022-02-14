from moha.system.operator.base import OperatorNames
from moha.system.hamiltonian.base import BaseHamiltonian
from moha.system.basis.slater_determinant import SlaterDeterminant
from moha.system.basis_set.ci_basis_set import CIBasisSet
import numpy as np
import copy
import itertools

class CIHamiltonian(object):
    """Configuration interaction Hamiltonian for a Schrodinger equation.

    Attributes
    ----------
    ham
        Chemical Hamiltonian.

    wfn
        Hartree Fock wavefunction.

    ci_basis_set
        Configuration interaction basis set.
    
    h1e : np.ndarray
        One electron integral.

    g2e : np.ndarray
        Two electron integral.

    Property
    -------
    npsin : int
        Number of spin orbitals

    Methods
    -------
    __init__(self,ham,wfn,ci_basis_set)
        Initialize the instance.
    
    assign_hamiltonian(self,ham)
        Assign the chemical Hamiltonian to the solver.
    
    assign_wavefunction(self,wfn)
        Assign the Hartree Fock wavefunction to the solver.

    generate_matrix(self,ci_basis_set)
        Generate CI Hamiltonian matrix.

    calculate_matrix_element(self,det1,det2)
        Calculate matrix element between two determinants.
    
    calculate_matrix_element_identical(self,det)
        Calculate matrix element whose configurations are identical.
    
    calculate_matrix_element_diff1(self,det1,det2):
        Calculate matrix element whose configuration differ by 1 spin orbitals.
    
    calculate_matrix_element_diff2(self,det1,det2):
        Calculate matrix element whose configuration differ by 2 spin orbitals.
    
    mixed_index_transformation(self,dic):
        Convert a configuration from dict form to list.
    """
    def __init__(self,ham,wfn,ci_basis_set):
        """Initialize the instance.

        Attributes
        ----------
        ham
            Chemical Hamiltonian.

        wfn
            Hartree Fock wavefunction.

        ci_basis_set
            Configuration interaction basis set.
        """        
        self.assign_hamiltonian(ham)
        self.assign_wavefunction(wfn)
        self.assign_ci_basis_set(ci_basis_set)
        self.assign_integral(ham,wfn)
    
    def assign_hamiltonian(self,ham):
        """Assign the chemical Hamiltonian to the solver.

        Parameters        
        ----------
        ham
            Chemical Hamiltonian.
        """
        self.ham = ham

    def assign_wavefunction(self,wfn):
        """Assign the Hartree Fock wavefunction to the solver.

        Parameters
        ----------
        wfn
            Hartree Fock wavefunction.
        """
        self.wfn = wfn

    def assign_ci_basis_set(self,ci_basis_set):
        """Assign the CI basis set to the instance.

        Parameters
        ----------
        ci_basis_set
            Configuration interaction basis set.
        """
        self.ci_basis_set = ci_basis_set
        
    def assign_integral(self,ham,wfn):
        """Assign the one and two electron integral in molecular spin orbital.
        
        Parameters
        ----------
        ham
            Chemical Hamiltonian.
        
        wfn
            Hartree Fock wavefunction.
        """
        C = wfn.coefficients
        
        Hcore = ham.operators[OperatorNames.Hcore].basis_transformation(C)
        h1e = Hcore.spin_orbital_basis_integral
        
        Eri = ham.operators[OperatorNames.Eri].basis_transformation(C)
        g2e = Eri.double_bar
        
        self.h1e = h1e
        self.g2e = g2e

    def integrate_wfn_sd(self, wfn, sd, wfn_deriv=None, ham_deriv=None):
        """Integrate the Hamiltonian with against a wavefunction and Slater determinant.
        """
        integral = 0.0
        print("{}=== {}".format(wfn.coefficients[0],self.calculate_matrix_element(sd,sd)))
        for i,deti in enumerate(wfn.basis_set.bases):
            print("{}=== {}".format(wfn.coefficients[i],self.calculate_matrix_element(sd,deti)))
            integral += wfn.coefficients[i]*self.calculate_matrix_element(sd,deti)
        return integral
    
    def generate_matrix(self,ci_basis_set):
        """Generate CI Hamiltonian matrix.

        Parameters
        ----------
        ci_basis_set : CIBasisSet
            Configuration interaction basis set.

        Raises
        ------
        TypeError
            If parameter ci_basis_set is not a CIBasisSet instance.

        Returns
        -------
        matrix : np.ndarray
            Configuration interaction Hamiltonian matrix.
        """
        if not isinstance(ci_basis_set,CIBasisSet):
            raise TypeError("Parameter ci_basis_set must be a CIBasisSet instance.")
        bases = ci_basis_set.bases
        num_det = ci_basis_set.size
        matrix = np.zeros((num_det,num_det))
        for i in range(num_det):
            for j in range(i+1):
                matrix[i,j] = self.calculate_matrix_element(bases[i],bases[j])
                matrix[j,i] = matrix[i,j]
        return matrix

    def calculate_matrix_element(self,det1,det2):
        """Calculate matrix element between two determinants.

        Parameters
        ----------
        det1 : SlaterDeterminant
            The first SlaterDeterminant instance.
        
        det2 : SlaterDeterminant
            The second SlaterDeterminant instance.
        
        Raises
        ------
        TypeError
            If det1 is not SlaterDeterminant instance.
            If det2 is not SlaterDeterminant instance.

        Returns
        -------
        value : float
            The integral value of CI Hamiltonian.
        """
        if not isinstance(det1,SlaterDeterminant):
            raise TypeError("Parameters det1 must be SlaterDeterminant instance.")
        if not isinstance(det2,SlaterDeterminant):
            raise TypeError("Parameters det2 must be SlaterDeterminant instance.")
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

        Parameters
        ----------
        det : SlaterDeterminant
            A SlaterDeterminant instance.
        
        Raises
        ------
        TypeError
            If det is not SlaterDeterminant instance.
        
        Returns
        -------
        value : float
            The integral value of CI Hamiltonian.
        """
        if not isinstance(det,SlaterDeterminant):
            raise TypeError("Parameters det must be SlaterDeterminant instance.")
        O1 = 0.0
        O2 = 0.0
        occ_mixed_index = self.mixed_index_transformation(det.occ_index)
        for i in occ_mixed_index:
            O1 += self.h1e[i,i]
        for i in occ_mixed_index:
            for j in occ_mixed_index:
                O2 += self.g2e[i,j,i,j]
        value = O1+0.5*O2
        return value

    def calculate_matrix_element_diff1(self,det1,det2):
        """Calculate matrix element whose configuration differ by 1 spin orbitals.

        Parameters
        ----------
        det1
            A SlaterDeterminant instance.
        
        det2
            A SlaterDeterminant instance.
        
        Raises
        ------
        TypeError
            If det1 is not SlaterDeterminant instance.
            If det2 is not SlaterDeterminant instance.
        
        Returns
        -------
        value : float
            The integral value of CI Hamiltonian.
        """
        if not isinstance(det1,SlaterDeterminant):
            raise TypeError("Parameters det1 must be SlaterDeterminant instance.")
        if not isinstance(det2,SlaterDeterminant):
            raise TypeError("Parameters det2 must be SlaterDeterminant instance.")
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
            O2 += self.g2e[a_list[0],k,c_list[0],k]
        value = O1+O2
        value *=sign
        return value

    def calculate_matrix_element_diff2(self,det1,det2):
        """Calculate matrix element whose configuration differ by 2 spin orbitals.

        Parameters
        ----------
        det1
            First SlaterDeterminant instance.
        
        det2
            Second SlaterDeterminant instance.
        
        Raises
        ------
        TypeError
            If det1 is not SlaterDeterminant instance.
            If det2 is not SlaterDeterminant instance.

        Returns
        -------
        value : float
            The integral value of CI Hamiltonian.
        """
        if not isinstance(det1,SlaterDeterminant):
            raise TypeError("Parameters det1 must be SlaterDeterminant instance.")
        if not isinstance(det2,SlaterDeterminant):
            raise TypeError("Parameters det2 must be SlaterDeterminant instance.")
        value = 0.0

        sign = det1.phase(det2)
        anni = det1.annihilation_list(det2)
        crea = det1.creation_list(det2)
        a_list = self.mixed_index_transformation(anni)
        c_list = self.mixed_index_transformation(crea)
        value += self.g2e[a_list[0],a_list[1],c_list[0],c_list[1]]
        value *=sign
        return value

    def mixed_index_transformation(self,dic):
        """Convert a configuration from dict form to list.

        Parameters
        ----------
        dic : dict
            A configuration in dict form.
        
        Returns
        -------
        l : list
            List representation of the configuration.
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
        
