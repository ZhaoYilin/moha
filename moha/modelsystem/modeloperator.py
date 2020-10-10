import numpy as np

class OneBodyTerm(object):
    """Operator for one body interaction"""

    def __init__(self,terms,cell1,site1,cell2,site2,parameter):
        """Initialize a OneBodyTerm instance
        
        Parameters
        ----------
        terms : list
            A list of length 2 for one body interaction

        cell1 : list
            Integer coordinate of first cell

        site1 : int
            Index of first site

        cell2 : list
            Integer coordinate of second cell

        site2 : int
            Index of second site
        
        parameter: float
            Parameter of the interaction

        """
        self.terms = terms
        self.cell1 = np.array(cell1)
        self.site1 = site1
        self.cell2 = np.array(cell2)
        self.site2 = site2
        self.parameter = parameter


class TwoBodyTerm(object):
    """Operator for two body interaction"""

    def __init__(self,terms,cell1,site1,cell2,site2,parameter):
        """Initialize a TwoBodyTerm instance
        
        Parameters
        ----------
        terms : list
            A list of two body interaction terms

        cell1 : list
            Integer coordinate of first and second cell

        site1 : int
            Index of first and second site

        cell2 : list
            Integer coordinate of third and fourth cell

        site2 : int
            Index of third and fourth site
        
        parameter: float
            Parameter of the interaction

        """
        self.terms = terms
        self.cell1 = np.array(cell1)
        self.site1 = site1
        self.cell2 = np.array(cell2)
        self.site2 = site2
        self.parameter = parameter

