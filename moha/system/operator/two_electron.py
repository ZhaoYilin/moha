from moha.system.operator.base import BaseOperator
import numpy as np

class TwoElectronOperator(BaseOperator):
    """
    """
    def __init__(self,name,nspatial,value=None):
        super(TwoElectronOperator,self).__init__(name,nspatial,value)

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

