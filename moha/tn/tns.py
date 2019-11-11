import numpy as np
import math
import copy
class TNS(object):
    """
    Base class for the general tensor network states(TNS): MPS, TTNW, PEPS, MERA, etc.

    TNS(U) := (V,E,H)
    V: set of vertices V= {U_i : i = 1,...,K}. K is the number of vertices.
    E: set of virtual indices dimension E = {m1,m2,...m_}. 
    H: set of physical indices H = {alpha_1,alpha_2,...,alpha_N}

    """
    def __init__(self,site,N,D):
        """
        Parameters:
            sites: object of sites class eg. SpinOneHalf SpinOne Spinless
            N: number of sites
            d: dimension of physical indices 
            D: maximum diemsion of virtual indices
        """
        super(TNS,self).__init__()
        self.N = N
        self.d = site.dim
        self.D = D

    def initializer(self):
        """
        initialize a tensor network states
        """

        raise NotImplementedError

class MPS(TNS):
    """
    Class for matrix product state(MPS)

    """

    def __init__(self,site,N,D,MPS_Type='MPS'):
        """
        Parameters:
            N: number of sites
            d: dimension of physical indices 
            D: maximum diemsion of virtual indices
        """
        self.d = site.dim
        self.N = N
        self.D = D
        self.MPS_Type = MPS_Type
        self._initializer(MPS_Type)
        
    def _initializer(self,MPS_Type):
        """
        Initialize a matrix product state

        A is a 3-index tensor, A[s,i,j]

           s
           |
        i -A- j
        
        [s] acts on the local Hilbert space
        [i,j] act on the virtual vonds

        Type: type of the matrix product state
            MPS: normal MPS
            uMPS: uniform MPS
            iMPS: infinite MPS

        BC: Boundary condition
            PBC periodic boundary condition
            OBC open boundary condition

        """


        N = self.N
        d = self.d
        D = self.D
        vertices = [None]*N

        if MPS_Type=='MPS':
            for i in range(N):
                vertices[i] = np.random.rand(d,min(D,d**i,d**(N-i)),min(D,d**(i+1),d**(N-i-1)))
            self.vertices = vertices

        elif MPS_Type=='uMPS':
            vertices[0] = np.random.rand(d,1,D)
            for i in range(1,N-1,1):
                vertices[i] = np.random.rand(d,D,D)
            vertices[N-1] = np.random.rand(d,D,1)
            self.vertices = vertices

        elif MPS_Type=='iMPS':
            for i in range(0,N,1):
                vertices[i] = np.random.rand(D,d,D)
            self.vertices = vertices

    def general_state_to_mps(self,tensor):
        """
        readin a general form of wavefunction and transfer it to MPS
        """
        N = self.N
        d = self.d
        D = self.D
        tensor = np.reshape(tensor,(1,)+tensor.shape+(1,)) 
        vertices = []
        
        bond_dim = [None]*(N+1)
        for i in range(N+1):
            bond_dim[i] = min(d**i,D,d**(N-i))
        for i in range(N-1):
            M = np.reshape(tensor,(np.prod(tensor.shape[0:2]),np.prod(tensor.shape[2:])))
            U,S,V = np.linalg.svd(M,full_matrices=0)
            S = np.diag(S)
            if len(S[0])>bond_dim[i+1]:
                U = U[:,0:bond_dim[i+1]]
                S = S[0:bond_dim[i+1]:,0:bond_dim[i+1]]
                V = V[0:bond_dim[i+1],:]
            U = np.reshape(U,(bond_dim[i],d,bond_dim[i+1]))
            V = np.reshape(V,[V.shape[0]]+[d]*int(math.log(V.shape[1],d))+[1])
            vertices.append(U)
            tensor = np.tensordot(S,V,axes=(1,0))
        vertices.append(tensor)
        for i in range(N):
            vertices[i] = np.transpose(vertices[i],(1,0,2))
        return vertices

    @property
    def general_state(self):
        """
        transfer wavefunction from MPS to general form
        """
        N = self.N
        d = self.d
        MPS = copy.deepcopy(self.vertices)
        for i in range(N):
            MPS[i] = np.transpose(MPS[i],(1,0,2))
        tensor = MPS[0]
        for i in range(1,N,1):
            tensor = np.tensordot(tensor,MPS[i],axes=(-1,0))
        tensor = np.trace(tensor,axis1=0,axis2=-1)

        return tensor

    def configuration(self,configuration):
        """
        Build a matrix product state for given configuration

        Parameters:
            configuration: a list of integer to represent a configuration
        
        """
        N = self.N
        d = self.d
        configuration_mps = copy.deepcopy(self.vertices)
        for i in range(N):
            for j in range(d):
                if j==configuration[i]:
                    pass
                else:
                    ml,md,mr = configuration_mps[i].shape
                    configuration_mps[i][:,j,:] = np.zeros(shape=(ml,mr))
        return configuration_mps

    def vertice_check(self,i):
        """
        Read-only access of vertice i

        """
        return self.vertices[i] 
    
    def vertice_update(self,A,i):
        """
        update the vertice i with new a tensor

        """
        if self.vertices[i].shape==A.shape:
            self.vertices[i] = copy.deepcopy(A)
        else:
            raise Exception('shape of vertice i and new tensor not match')

    def vertice_swap(self,i,j):
        """
        Swap vertice i and j

        """
        if self.vertices[i].shape==self.vertices[j].shape:
            self.vertices[i],self.vertices[j] = copy.deepcopy(self.vertices[j]),copy.deepcopy(self.vertices[i])
        else:
            raise Exception('shape of vertice i and j not match')


    def guage(self,i,j):
        """
        MPS are not quique, we can defined up to gauge on the virtual bond

            Parameters:
                i: index of ith vertice
                j: index of jth vertice
        """
        if abs(i-j)!=1:
            raise ValueError('i and j must be adjacent')
        #todo add more restrict condition
        D = self.vertices[i].shape[2]
        G = np.random.rand(D,D)
        G_inv = np.linalg.inv(G)
        self.vertices[i] = np.tensordot(self.vertices[i],G,axes=(2,0))
        self.vertices[j] = np.tensordot(G_inv,self.vertices[j],axes=(1,0))

    @property
    def norm(self):
        """
        norm of a MPS, return back a scalar
        """
        pass
    def norm_vertice(self,i,direction):
        """
        Get the norm of vertice i

        """
        pass
    def norm_left(self,index):
        """
        Calculate the norm of the MPS on the left side of vertice i 
        """
        pass
    def norm_right(self,index):
        """
        Calculate the norm of the MPS on the right side of vertice i 
        """
        pass

    def vertice_decomposition(self,index,direction):
        """
        index: index of the vertice
        """
        m_l,m_d,m_r = self.vertices[index].shape

        if direction=='left':
            A = copy.deepcopy(self.vertices[index])
            M = np.reshape(A,(m_l*m_d,m_r))
            U,S,V = np.linalg.svd(M,full_matrices=0)
            return U,S,V

        elif direction=='right':
            A = copy.deepcopy(self.vertices[index])
            M = np.reshape(A,(m_l,m_d*m_r))
            U,S,V = np.linalg.svd(M,full_matrices=0)
            return U,S,V

        else:
            raise ValueError('direction must be either left or right')

class PEPS(TNS):
    pass

class MERA(TNS):
    pass
