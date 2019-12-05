import numpy as np
from moha.tools.orbital_order.wblossom import find_maximum_matching


class BOTTOMUP:
    def __init__(self,adj):
        self.adj = adj
    
    def adj_convert(self,adj):
        """Convert adj matrix to list format

        """
        adj_list = []
        for i in range(adj.shape[0]):
            for j in range(adj.shape[0]):
                if i<j:
                    l = tuple([i]+[j]+[adj[i][j]])
                    adj_list.append(l)
        return adj_list

    def match(self,adj_list):
        """
        """

        match = find_maximum_matching(adj_list)
        return match

    def adj_update(self,adj_old):
        d = adj_old.shape[0]
        adj_new = np.zeros(adj_old.shape)

        for i in range(d):
            for j in range(d):
                print(d)
                adj_new[i][j] = adj_old[2*i-1][2*j-1]+adj_old[2*i-1][2*j]+adj_old[2*i][2*j-1]+adj_old[2*i][2*j]
       
        return adj_new 
