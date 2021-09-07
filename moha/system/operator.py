from abc import ABC,abstractmethod,abstractproperty
import numpy as np

class BaseOperator(ABC):
    """Base Operator

    Attributes
    ----------
    name : str
        Name of the operator

    nspatial : int
        Number of the spatial orbitals

    value : np.array
        Value of the operator

    Abstract Property
    -----------------
    dim : int
        Dimension of the operator

    npsin : int
        Number of spin orbitals
    
    Methods
    -------
    __init__(self,name,nspatial,value)
        Initialize the Operator.        
    
    Abstract Method
    ---------------
    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.
    
    assign_value(self,value)
        Assign integral value for the operator.
    """

    def __init__(self,name,nspatial,value=None):
        """Initialize the base operator

        Parameters
        ----------
        name : str
            Name of the operator.

        nspatial : int
            Number of the spatial orbitals for the operator

        value : np.array
            Integral value for the operator.

        """
        self.assign_name(name)
        self.assign_nspatial(nspatial)
        self.assign_value(value)
        
    @abstractproperty
    def dim(self):
        """Return the dimension of operator.

        Returns
        -------
        dim : int
            Dimension of operator.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError

    @abstractproperty
    def nspin(self):
        """Return the number of spin orbitals.

        Returns
        -------
        nspin : int
            Number of spin orbitals.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError

    @abstractmethod
    def assign_name(self,name):
        """Assign name to the operator.

        Parameters
        ----------
        name : str
            Name of the operator.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError 

    @abstractmethod
    def assign_nspatial(self,nspatial):
        """Assign number of spatial orbitals for the operator.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbitals for the operator.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError 
    
    @abstractmethod
    def assign_value(self,value):
        """Assign integral value for the operator.

        Parameters
        ----------
        value : np.array
            Integral value for the operator.

        Raises
        ------
        NotImplementedError
            If called.
        """
        raise NotImplementedError 

class ZeroElectronOperator(BaseOperator):
    """Base Operator

    Attributes
    ----------
    name : str
        Name of the operator

    nspatial : int
        Number of the spatial orbitals

    value : np.array
        Value of the operator

    Property
    --------
    dim : int
        Dimension of the operator

    npsin : int
        Number of spin orbitals
    
    Methods
    -------
    __init__(self,name,nspatial,value)
        Initialize the Operator.        
    
    build(cls,name,molecule)
        Build the operator instance
    
    assign_name(self,name)
        Assign name to the operator.

    assign_nspatial(self,nspatial)
        Assign number of spatial orbitals for the operator.
    
    assign_value(self,value)
        Assign integral value for the operator.
    """
    def __init__(self,name,nspatial,value):
        self.assign_name(name)
        self.assign_nspatial(nspatial)
        self.assign_value(value)

    @classmethod
    def build(cls,name,molecule):
        value = 0.
        if name=='nuclear_repulsion':
            for i,A in enumerate(molecule.atoms):
                for j,B in enumerate(molecule.atoms):
                    if j>i:
                        value += (A.number*B.number)/molecule.bond_length(i,j)
        else:
            raise ValueError('Unknown operator or operator not implemented')
        operator = cls(name,nspatial,value)
        return operator

    @property        
    def dim(self):
        """Return the dimension of operator.

        Returns
        -------
        dim : int
            Dimension of operator.

        Raises
        ------
        TypeError
            If integral value is not None or np.ndarray.        
        """
        if isinstance(self.value,np.ndarray):
            self.dim = self.value.shape[0]
        elif isinstance(self.value,None):
            self.dim = self.nspatial
        else:
            raise TypeError("value must be None or np.ndarray")


    @property        
    def nspin(self):
        """Return the number of spin orbitals.

        Returns
        -------
        nspin : int
            Number of spin orbitals.
        """
        self.nspin = nspatial*2

    def assign_name(self,name):
        """Assign name to the operator.

        Parameters
        ----------
        name : str
            Name of the operator.

        Raises
        ------
        TypeError
            If name of operator is not str.
        """
        if not isinstance(name, str):
            raise TypeError("Name of oeprator must be str")
        self.name = name

    def assign_nspatial(self,nspatial):
        """Assign number of spatial orbitals for the operator.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbitals for the operator.

        Raises
        ------
        TypeError
            If number of electrons is not an integer.
        ValueError
            If number of electrons is not a positive number.
        """
        if not isinstance(nspatial, int):
            raise TypeError("Number of spatial orbitals must be an integer")
        if nspatial <= 0:
            raise ValueError("Number of spatial orbitals must be a positive integer")
        self.nspatial = nspatial
    
    def assign_value(self,value):
        """Assign integral value for the operator.

        Parameters
        ----------
        value : np.array
            Integral value for the operator.

        Raises
        ------
        TypeError
            If integral value is not np.ndarray.
        """
        if not isinstance(value, np.ndarray):
            raise TypeError("Integral value must be np.ndarray")
        self.value = value

class OneElectronOperator(BaseOperator):
    """
    """
    def __init__(self,name,dim,value):
        super(OneElectronOperator,self).__init__(name,dim,value)

    @classmethod
    def build(cls,name,molecule,basis_set,*args):
        Norb = basis_set.size
        value = np.zeros((Norb,Norb))     
        if name=='overlap':
            from .integral.overlap import overlap
            for i,oi in enumerate(basis_set.basis):
                for j,oj in enumerate(basis_set.basis):
                    if i>=j:
                        value[i][j] = overlap(oi,oj)
        elif name == 'kinetic':
            from .integral.kinetic import kinetic
            for i,oi in enumerate(basis_set.basis):
                for j,oj in enumerate(basis_set.basis):
                    if i>=j:
                        value[i][j] = kinetic(oi,oj)
        elif name == 'nuclear_attraction':
            from .integral.nuclear_attraction import nuclear_attraction
            for i,oi in enumerate(basis_set.basis):
                for j,oj in enumerate(basis_set.basis):
                    if i>=j:
                        for atom in molecule.atoms:
                            value[i][j] += -1*atom.number*nuclear_attraction(oi,oj,atom.coordinate)
        elif name == 'multipole_moment':
            from .integral.multipole_moment import multipole_moment
            for i,oi in enumerate(basis_set.basis):
                for j,oj in enumerate(basis_set.basis):
                    if i>=j:
                        value[i][j] = multipole_moment(oi,oj,args[0],args[1])
        elif name == 'differential':
            from .integral.differential import differential
            for i,oi in enumerate(basis_set.basis):
                for j,oj in enumerate(basis_set.basis):
                    if i>=j:
                        value[i][j] = differential(oi,oj,args[0])
        elif name == 'linear_momentum':
            from .integral.linear_momentum import linear_momentum
            value = value.astype(np.complex)
            for i,oi in enumerate(basis_set.basis):
                for j,oj in enumerate(basis_set.basis):
                    if i>=j:
                        value[i][j] = linear_momentum(oi,oj)
        elif name == 'angular_momentum':
            from .integral.angular_momentum import angular_momentum
            value = value.astype(np.complex)
            for i,oi in enumerate(basis_set.basis):
                for j,oj in enumerate(basis_set.basis):
                    if i>=j:
                        value[i][j] = angular_momentum(oi,oj)
        else:
            raise ValueError('Unknown operator or operator not implemented')

        value = value + value.T - np.diag(value.diagonal())
        operator = cls(name,basis_set.size,value)
        return operator

    @property        
    def dim(self):
        """Return the dimension of operator.

        Returns
        -------
        dim : int
            Dimension of operator.

        Raises
        ------
        TypeError
            If integral value is not None or np.ndarray.        
        """
        if isinstance(self.value,np.ndarray):
            self.dim = self.value.shape[0]
        elif isinstance(self.value,None):
            self.dim = self.nspatial
        else:
            raise TypeError("value must be None or np.ndarray")


    @property        
    def nspin(self):
        """Return the number of spin orbitals.

        Returns
        -------
        nspin : int
            Number of spin orbitals.
        """
        self.nspin = nspatial*2

    def assign_name(self,name):
        """Assign name to the operator.

        Parameters
        ----------
        name : str
            Name of the operator.

        Raises
        ------
        TypeError
            If name of operator is not str.
        """
        if not isinstance(name, str):
            raise TypeError("Name of oeprator must be str")
        self.name = name

    def assign_nspatial(self,nspatial):
        """Assign number of spatial orbitals for the operator.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbitals for the operator.

        Raises
        ------
        TypeError
            If number of electrons is not an integer.
        ValueError
            If number of electrons is not a positive number.
        """
        if not isinstance(nspatial, int):
            raise TypeError("Number of spatial orbitals must be an integer")
        if nspatial <= 0:
            raise ValueError("Number of spatial orbitals must be a positive integer")
        self.nspatial = nspatial
    
    def assign_value(self,value):
        """Assign integral value for the operator.

        Parameters
        ----------
        value : np.array
            Integral value for the operator.

        Raises
        ------
        TypeError
            If integral value is not np.ndarray.
        """
        if not isinstance(value, np.ndarray):
            raise TypeError("Integral value must be np.ndarray")
        self.value = value

    def basis_transformation(self,C):
        """
        C:coefficent matrix
        """
        if type(C) is np.ndarray:
            value_prime = np.zeros((self.dim,self.dim))
            for i in range(self.dim): 
                for j in range(self.dim): 
                    for k in range(self.dim):
                        for l in range(self.dim):
                            value_prime[i][j] += C[k][j]*C[l][i]*self.value[k][l]
            self.value = value_prime
            return value_prime

        elif type(C) is dict:
            value_prime = {'alpha':np.zeros((self.dim,self.dim)),'beta':np.zeros((self.dim,self.dim))}
            for spin in ['alpha','beta']:
                for i in range(self.dim): 
                    for j in range(self.dim): 
                        for k in range(self.dim):
                            for l in range(self.dim):
                                value_prime[spin][i,j] += C[spin][k,j]*C[spin][l,i]*self.value[spin][k,l]
            self.value = value_prime
            return value_prime

    @property
    def spin(self):
        if type(self.value) is np.ndarray:
            value_prime = np.zeros((self.dim*2,self.dim*2))
            for i in range(self.dim*2): 
                for j in range(self.dim*2): 
                    if i%2==j%2:
                        value_prime[i,j]  = self.value[int(i/2),int(j/2)]
            return value_prime

        elif type(self.value) is dict:
            value_prime = np.zeros((self.dim*2,self.dim*2))
            for i in range(self.dim*2): 
                for j in range(self.dim*2): 
                    if i%2==j%2==0:
                        value_prime[i,j]  = self.value['alpha'][int(i/2),int(j/2)]
                    elif i%2==j%2==1:
                        value_prime[i,j]  = self.value['beta'][int(i/2),int(j/2)]
            return value_prime



class TwoElectronOperator(BaseOperator):
    """
    """
    def __init__(self,name,dim,value):
        super(TwoElectronOperator,self).__init__(name,dim,value)

    @classmethod
    def build(cls,name,molecule,basis_set):
        if name=='electron_repulsion':
            from .integral.electron_repulsion import electron_repulsion
            value = []     
            for i,oi in enumerate(basis_set.basis):
                for j,oj in enumerate(basis_set.basis[0:i+1]):
                    for k,ok in enumerate(basis_set.basis[0:i+1]):
                        for l,ol in enumerate(basis_set.basis[0:k+1]):
                            ij = i*(i+1)*0.5+j
                            kl = k*(k+1)*0.5+l
                            if ij>=kl:
                                value.append(electron_repulsion(oi,oj,ok,ol))
        else:
            raise ValueError('Unknown operator or operator not implemented')
        value = np.array(value)
        operator = cls(name,basis_set.size,value)
        return operator

    @property        
    def dim(self):
        """Return the dimension of operator.

        Returns
        -------
        dim : int
            Dimension of operator.

        Raises
        ------
        TypeError
            If integral value is not None or np.ndarray.        
        """
        if isinstance(self.value,np.ndarray):
            self.dim = self.value.shape[0]
        elif isinstance(self.value,None):
            self.dim = self.nspatial
        else:
            raise TypeError("value must be None or np.ndarray")


    @property        
    def nspin(self):
        """Return the number of spin orbitals.

        Returns
        -------
        nspin : int
            Number of spin orbitals.
        """
        self.nspin = nspatial*2

    def assign_name(self,name):
        """Assign name to the operator.

        Parameters
        ----------
        name : str
            Name of the operator.

        Raises
        ------
        TypeError
            If name of operator is not str.
        """
        if not isinstance(name, str):
            raise TypeError("Name of oeprator must be str")
        self.name = name

    def assign_nspatial(self,nspatial):
        """Assign number of spatial orbitals for the operator.

        Parameters
        ----------
        nspatial : int
            Number of spatial orbitals for the operator.

        Raises
        ------
        TypeError
            If number of electrons is not an integer.
        ValueError
            If number of electrons is not a positive number.
        """
        if not isinstance(nspatial, int):
            raise TypeError("Number of spatial orbitals must be an integer")
        if nspatial <= 0:
            raise ValueError("Number of spatial orbitals must be a positive integer")
        self.nspatial = nspatial
    
    def assign_value(self,value):
        """Assign integral value for the operator.

        Parameters
        ----------
        value : np.array
            Integral value for the operator.

        Raises
        ------
        TypeError
            If integral value is not np.ndarray.
        """
        if not isinstance(value, np.ndarray):
            raise TypeError("Integral value must be np.ndarray")
        self.value = value

    def tensor_format(self):
        """
        transform electron repulsion integral from list format to tensor format
        """
        from .auxiliary import eint
        Norb = self.dim
        Eri_tensor = np.zeros((Norb,Norb,Norb,Norb))
        for i in range(Norb):
            for j in range(Norb):
                for k in range(Norb):
                    for l in range(Norb):
                        Eri_tensor[i,j,k,l] = self.value[eint(i,j,k,l)]
 
        return Eri_tensor

    def basis_transformation(self,C):
        """
        transform electron repulsion integral from atomic orbtial to molecular orbtial
        C:coefficent matrix
        """
        from .auxiliary import eint
        Norb = self.dim
        if type(C) is np.ndarray:
            Eri_tensor = np.zeros((Norb,Norb,Norb,Norb))
            for i in range(Norb):
                for j in range(Norb):
                    for k in range(Norb):
                        for l in range(Norb):
                            Eri_tensor[i,j,k,l] = self.value[eint(i,j,k,l)]
 
            temp = np.zeros((Norb,Norb,Norb,Norb))
            temp2 = np.zeros((Norb,Norb,Norb,Norb))
            temp3= np.zeros((Norb,Norb,Norb,Norb))
            Eri = np.zeros((Norb,Norb,Norb,Norb))
            for i in range(Norb):
                for m in range(Norb):
                    temp[i,:,:,:] += C[m,i]*Eri_tensor[m,:,:,:]
                for j in range(Norb):
                    for n in range(Norb):
                        temp2[i,j,:,:] += C[n,j]*temp[i,n,:,:]
                    for k in range(Norb):
                        for o in range(Norb):
                            temp3[i,j,k,:] += C[o,k]*temp2[i,j,o,:]
                        for l in range(Norb):
                            for p in range(Norb):
                                Eri[i,j,k,l] += C[p,l]*temp3[i,j,k,p]
            self.value = Eri
            return Eri

        elif type(C) is dict:
            raise ValueError('basis type should be spatial')
    @property
    def coulomb(self):
        value_prime = np.zeros((self.dim,self.dim))
        for i in range(self.dim):
            for j in range(self.dim):
                value_prime[i,j] = self.value[i,i,j,j]
        return value_prime

    @property
    def exchange(self):
        value_prime = np.zeros((self.dim,self.dim))
        for i in range(self.dim):
            for j in range(self.dim):
                value_prime[i,j] = self.value[i,j,i,j]
        return value_prime

    @property
    def spin(self):
        if type(self.value) is np.ndarray:
            value_prime = np.zeros((self.dim*2,self.dim*2,self.dim*2,self.dim*2))
            for i in range(self.dim*2): 
                for j in range(self.dim*2): 
                    for k in range(self.dim*2): 
                        for l in range(self.dim*2): 
                            if i%2==k%2 and j%2==l%2:
                                value_prime[i,j,k,l]  = self.value[int(i/2),int(k/2),int(j/2),int(l/2)]
            return value_prime

        elif type(self.value) is dict:
            raise ValueError('basis type should be spatial')

    @property
    def double_bar(self):
        """
        antisymmetrized two electron integral
        """
        value = self.spin
        if type(self.value) is np.ndarray:
            value_prime = np.zeros((self.dim*2,self.dim*2,self.dim*2,self.dim*2))
            for i in range(self.dim*2): 
                for j in range(self.dim*2): 
                    for k in range(self.dim*2): 
                        for l in range(self.dim*2): 
                            value_prime[i,j,k,l]  = value[i,j,k,l] - value[i,j,l,k]
            return value_prime

        elif type(self.value) is dict:
            raise ValueError('basis type should be spatial')

