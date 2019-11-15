import numpy as np

class Operator(object):
    """
    """
    def __init__(self,name,dim,value):
        self.name = name
        self.dim = dim
        self.value = value

class CreationOperator(object):
    def __init__(self):
        pass

class AnnihilationOperator(object):
    def __init__(self):
        pass

class ZeroElectronOperator(Operator):
    """
    """
    def __init__(self,name,value):
        self.name = name
        self.value = value

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
        operator = cls(name,value)
        return operator
        
class OneElectronOperator(Operator):
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

    def spin(self):
        if type(self.value) is np.ndarray:
            value_prime = np.zeros((self.dim*2,self.dim*2))
            for i in range(self.dim*2): 
                for j in range(self.dim*2): 
                    if i%2==j%2:
                        value_prime[i,j]  = self.value[i/2,j/2]
            return value_prime

        elif type(self.value) is dict:
            value_prime = np.zeros((self.dim*2,self.dim*2))
            for i in range(self.dim*2): 
                for j in range(self.dim*2): 
                    if i%2==j%2==0:
                        value_prime[i,j]  = self.value['alpha'][i/2,j/2]
                    elif i%2==j%2==1:
                        value_prime[i,j]  = self.value['beta'][i/2,j/2]
            return value_prime



class TwoElectronOperator(Operator):
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

class Hamiltonian(object):
    """
    """
    def __init__(self, dim, operators={}):
        self.dim = dim
        self.operators = operators

    def add_operator_by_value(self,operator_name,type,value):
        if str(operator_name) in self.operators: 
            raise Exception("Operator name exists already")
        else:
            if type==0:
                # warning
                self.operators[str(operator_name)] = ZeroElectronOperator(operator_name,value)
            elif type==1:
                self.operators[str(operator_name)] = OneElectronOperator(operator_name,self.dim,value)
            elif type==2:
                self.operators[str(operator_name)] = TwoElectronOperator(operator_name,self.dim,value)
            
            else:
                raise ValueError('operators can only be zero one or two electrons operator')

    def add_operator(self,molecule,basis_set,operator_name,type):
        """Adds an operator to the hamiltonian.
    
        Parameters
        ----------
        operator_name : string
        The operator name.
        type : int
            type of operator 
        """
        if str(operator_name) in self.operators:
    	    raise ValueError("Operator name exists already")
        else:
            if type==0:
                self.operators[str(operator_name)] = ZeroElectronOperator.build(operator_name,molecule)
            elif type==1:
                self.operators[str(operator_name)] = OneElectronOperator.build(operator_name,molecule,basis_set)
            elif type==2:
                self.operators[str(operator_name)] = TwoElectronOperator.build(operator_name,molecule,basis_set)
            else:
                raise ValueError('operators can only be zero one or two electrons operator')


    @classmethod
    def build(cls,molecule,basis_set,terms={'nuclear_repulsion':0,'overlap':1,'kinetic':1,'nuclear_attraction':1,'electron_repulsion':2}):
        dim = basis_set.size
        hamiltonian = cls(dim)
        for operator_name,type in terms.items():
            hamiltonian.add_operator(molecule,basis_set,operator_name,type)
        return hamiltonian

    def basis_transfer(self,C,operator_name):
        if str(operator_name) not in self.operators:
            raise ValueError("Operator name not exists already")
        else:
            return self.operators[operator_name].basis_transfer(C)

