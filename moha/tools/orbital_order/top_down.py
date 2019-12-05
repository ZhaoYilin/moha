import numpy as np
from numpy import linalg as LA
from statistics import median

class TOPDOWN:

    def __init__(self, adjacency,type='zero'):
        """Divide vector into two sets. 

        Parameters
        ----------
        adjacency: numpy array 
            The adjacency matrix
        type: str
            type of division, it can be 'zero', 'median' or 'gap'.

        """

        self.adj = adjacency
        self.type = type
        self.dim = adjacency.shape[0]
        self.list = list(range(self.dim))
        self.order = []
        self.lap = self.laplacian(adjacency)


    def laplacian(self,adj):
        dim = self.dim
        lap = np.zeros((dim,dim))
        for i in range(dim):
            for j in range(dim):
                if i==j:
                    lap[i][j] = adj[i][j]
                else:
                    lap[i][j] = -adj[i][j]
        return lap
    
    def fiedler_vector(self, matrix):
        """The eigenvector corresponding to the second smallest eigenvalue
        of the Laplacian matrix.
        
        Parameters
        ----------
        matrix : numpy array
            The Laplacian matrix of a graph.
        """
        eigenValues, eigenVectors = LA.eigh(matrix)
        return eigenVectors[1]


    def divide(self, lap, abs_order):
        """Divide vector into two sets. 

        Parameters
        ----------
        lap: numpy array 
            The Laplacian matrix
        abs_order: list
            List of order using the index from original matrix.
            
        """
        fv = self.fiedler_vector(lap)
        if self.type=='zero':
            m=0
        elif self.type=='median':
            m = median(fv)
        elif self.type=='gap':
            m = np.argmax(np.diff(sorted(fv)))
        else:
            raise Exception('type of dividsion can be only zero median or gap')

        S1 = [i for i, x in enumerate(fv) if x < m]
        S2 = [i for i, x in enumerate(fv) if x >= m]
        S1 = [abs_order[i] for i in S1]
        S2 = [abs_order[i] for i in S2]
                
        return S1,S2
    

    def sub_block(self, ls):
        """
        """
        return self.lap[np.ix_(ls,ls)]


    def core(self,lap,abs_order):
        """
        """
        if len(abs_order)==1:
            self.order += abs_order
        else:
            S1,S2 = self.divide(lap,abs_order)
            A1,A2 = self.sub_block(S1), self.sub_block(S2)

            return self.core(A1,S1),self.core(A2,S2)


    def run(self):
        self.core(self.lap,self.list)
        return self.order












