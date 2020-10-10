import numpy as np
import copy

class SlaterDeterminant(object):
    """Slater Determinant Basis.

    Attributes
    ----------
    configuration : dict
        Binary representation of the determinants.
        {'alpha':[1,1,0,0],'beta':[1,1,0,0]}

    Properties
    ----------
    configuration_alpha : list
        configuraiton of alpha spin orbitals.
    
    configuration_beta : list
        configuraiton of beta spin orbitals.
    
    occ_index : list
        occupation index of alpha and beta spin orbitals.
    
    Methods
    -------
    __init__(self,configuration={})
        Initialize the instance.
    
    assign_configuration(self,configuration)
        Assign configuration to the Slater determinant.
    
    degree_of_excitation_alpha(self,det)
        Compute the degree of excitation between D1 and D2 for alpha spin.
    
    degree_of_excitation_beta(self,det)
        Compute the degree of excitation between D1 and D2 for beta spin.
    
    degree_of_excitation(self,det)
        Compute the degree of excitation between D1 and D2.
    
    annihilation_list_alpha(self,det)
        Identifying the annihilated alpha spin orbitals.
    
    annihilation_list_beta(self,det)
        Identifying the annihilated beta spin orbitals.
    
    creation_list_alpha(self,det)
        Identifying the created alpha spin orbitals.
    
    creation_list_beta(self,det)
        Identifying the created beta spin orbitals.
    
    annihilation_list(self,det)
        Identifying the annihilated spin orbitals.
    
    creation_list(self,det)
        Identifying the created spin orbitals.
    
    phase(self,det)
        The phase is calculated as −1^{Nperm}.This number is equal
        to the number of occupied spin-orbitals between these two positions.
    """

    def __init__(self,configuration={}):
        """Initialize the instance.

        Parameters
        ----------
        configuration : dict
            Binary representation of the Slater determinant.
        """
        self.assign_configuration(configuration)

    @property
    def configuration_alpha(self):
        """Configuraiton of alpha spin orbitals.
        
        Returns
        -------
        configuration_alpha : list
            configuraiton of alpha spin orbitals.
        """
        return self.configuration['alpha']

    @property
    def configuration_beta(self):
        """Configuraiton of beta spin orbitals.
        
        Returns
        -------
        configuration_beta : list
            configuraiton of beta spin orbitals.
        """
        return self.configuration['beta']

    @property
    def occ_index(self):
        """Occupation index of alpha and beta spin orbitals.
        
        Returns
        -------
        occ_index : dict
            Occupation index of alpha and beta spin orbitals.
        """
        occ_alpha = []
        occ_beta = []
        for index,item in enumerate(self.configuration['alpha']):
            if item==1:
                occ_alpha.append(index)
        for index,item in enumerate(self.configuration['beta']):
            if item==1:
                occ_beta.append(index)
        return {"alpha":occ_alpha,"beta":occ_beta}
    
    def assign_configuration(self,configuration):
        """Assign configuration to the Slater determinant.
        
        Parameters
        ----------
        configuration : dict
            Binary representation of the Slater determinant.

        Raises
        ------
        TypeError
            If configuration is not a dictionary.
        """
        if not isinstance(configuration,dict):
            raise TypeError("Configuration must be a dictionary.")
        self.configuration = configuration
    
    def degree_of_excitation_alpha(self,det):
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
        if not isinstance(det,SlaterDeterminant):
            raise TypeError("Parameter det must be a SlaterDeterminant instance.")
        diff = list(np.array(det.configuration['alpha']) - np.array(self.configuration['alpha']))
        degree = int(sum(map(abs,diff))/2)
        return degree
    
    def degree_of_excitation_beta(self,det):
        """Compute the degree of excitation between D1 and D2 for beta spin.

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
            Degree of beta spin excitation between two Slater determinant.
        """
        if not isinstance(det,SlaterDeterminant):
            raise TypeError("Parameter det must be a SlaterDeterminant instance.")
        diff = list(np.array(det.configuration['beta']) - np.array(self.configuration['beta']))
        degree = int(sum(map(abs,diff))/2)
        return degree

    def degree_of_excitation(self,det):
        """Compute the degree of excitation between D1 and D2.

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
            Degree of excitation between two Slater determinant.
        """
        if not isinstance(det,SlaterDeterminant):
            raise TypeError("Parameter det must be a SlaterDeterminant instance.")
        alpha_degree = self.degree_of_excitation_alpha(det)
        beta_degree = self.degree_of_excitation_beta(det)
        degree = alpha_degree + beta_degree
        return degree

    def annihilation_list_alpha(self,det):
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
        if not isinstance(det,SlaterDeterminant):
            raise TypeError("Parameter det must be a SlaterDeterminant instance.")
        diff = np.array(det.configuration['alpha']) - np.array(self.configuration['alpha'])
        if np.sum(np.abs(diff)) == 0: 
            index = []
        else:
            index = np.where(diff==-1)[0].tolist()
        return index


    def annihilation_list_beta(self,det):
        """Identifying the annihilated beta spin orbitals.

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
            Index of annihilation beta spin.
        """
        if not isinstance(det,SlaterDeterminant):
            raise TypeError("Parameter det must be a SlaterDeterminant instance.")
        diff = np.array(det.configuration['beta']) - np.array(self.configuration['beta'])
        if np.sum(np.abs(diff)) == 0: 
            index = []
        else:
            index = np.where(diff==-1)[0].tolist()
        return index

    def creation_list_alpha(self,det):
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
        if not isinstance(det,SlaterDeterminant):
            raise TypeError("Parameter det must be a SlaterDeterminant instance.")
        diff = np.array(det.configuration['alpha']) - np.array(self.configuration['alpha'])
        if np.sum(np.abs(diff)) == 0: 
            index = []
        else:
            index = np.where(diff==1)[0].tolist()
        return index

    def creation_list_beta(self,det):
        """Identifying the created beta spin orbitals.
        
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
            Index of creation beta spin.
        """
        if not isinstance(det,SlaterDeterminant):
            raise TypeError("Parameter det must be a SlaterDeterminant instance.")
        diff = np.array(det.configuration['beta']) - np.array(self.configuration['beta'])
        if np.sum(np.abs(diff)) == 0: 
            index = []
        else:
            index = np.where(diff==1)[0].tolist()
        return index

    def annihilation_list(self,det):
        """Identifying the annihilated spin orbitals.
        
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
            Index of annihilation spin.
        """
        if not isinstance(det,SlaterDeterminant):
            raise TypeError("Parameter det must be a SlaterDeterminant instance.")
        alpha_list = self.annihilation_list_alpha(det)
        beta_list = self.annihilation_list_beta(det)
        index = {'alpha':alpha_list, 'beta':beta_list}
        return index
    
    def creation_list(self,det):
        """Identifying the created spin orbitals.
        
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
            Index of creation spin.
        """
        if not isinstance(det,SlaterDeterminant):
            raise TypeError("Parameter det must be a SlaterDeterminant instance.")
        alpha_list = self.creation_list_alpha(det)
        beta_list = self.creation_list_beta(det)
        index = {'alpha':alpha_list, 'beta':beta_list}
        return index

    def phase(self,det):
        """The phase is calculated as −1^{Nperm}.This number is equal
        to the number of occupied spin-orbitals between these two positions.
        
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
        sign : int
            Phase for the evaluation between two Slater determinant.
        """
        if not isinstance(det,SlaterDeterminant):
            raise TypeError("Parameter det must be a SlaterDeterminant instance.")
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

