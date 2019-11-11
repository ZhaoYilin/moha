import numpy as np

class Site(object):
    """A class for general single site
    
    Use this class to create a single site object. The site comes with identity
    operator for a given dimension. To build specific site, additional operators
    need be add with add_operator method.
    
    """
    def __init__(self, dim):
        """Creates an empty site of dimension dim.

        Parameters
        ----------

        dim : an int 
            Size of the Hilbert space for single site. The dimension must be at 
            least 1. A site of dim = 1 is trival which represents the vaccum

        operators : a dictionary of string and numpy array (with ndim = 2).
            Operators for the site.

        """
        super(Site, self).__init__()
        self.dim = dim
        self.states = {}
        self.operators = { "id" : np.eye(self.dim, self.dim) }
    
    def add_operator(self, operator_name):
        """Adds an operator to the site with zero matrix.
    
        Parameters
        ----------
        operator_name : string
            The operator name.
	
        """
        self.operators[str(operator_name)] = np.zeros((self.dim, self.dim))

    def add_state(self, state_name):
        """Adds an state to the site with zero list.
    
        Parameters
        ----------
        operator_name : string
            The operator name.
	
        """
        self.states[str(state_name)] = np.zeros(self.dim)

class SpinlessFermionSite(Site):
    """A site for spinless fermion models.
    
    Use this class for spinless fermion sites. The Hilbert space is ordered
    such as:
    
    - the first state is empty site
    - the second state is occupied site.

    Notes
    -----
    Postcondition : The site has already built-in the operators for
                    c, c_dag, n.

    """
    def __init__(self):
        """Creates the spin one-half site.

        Notes
        -----
        Postcond : the dimension is set to 2

        """
        super(SpinlessFermionSite, self).__init__(2)
        # add the operators
        self.add_operator("c")
        self.add_operator("c_dag")
        self.add_operator("n")
        # for clarity
        c = self.operators["c"]
        c_dag = self.operators["c_dag"]
        n = self.operators["n"]
        # set the matrix elements different from zero to the right values
        c[0, 1] = 1
        c_dag[1, 0] = 1
        n[1, 1] = 1
        # add the states
        self.add_state("empty")
        self.add_state("occupied")
        # for clarity
        state_empty = self.states["empty"]
        state_occupied = self.states["occupied"]
        # set the list elements different from zero to the right values
        state_empty[0] = 1.0
        state_occupied[1] = 1.0


class SpinOneHalfSite(Site):
    """A site for spin 1/2 models.
    
    Use this class for spin one-half sites. The Hilbert space is ordered
    such as the first state is the spin down, and the second state is the
    spin up. 

    Notes
    -----
    Postcondition : The site has already built-in the spin operators for
                    s_x, s_y, s_z, s_p, s_m.

    """
    def __init__(self):
        """Creates the spin one-half site.

        Notes
        -----
        Postcond : the dimension is set to 2

        """
        super(SpinOneHalfSite, self).__init__(2)
        # add the operators
        self.add_operator("s_x")
        self.add_operator("s_y")
        self.add_operator("s_z")
        self.add_operator("s_p")
        self.add_operator("s_m")
        # for clarity
        s_x = self.operators["s_x"]
        s_y = self.operators["s_y"]
        s_y = s_y.astype(np.complex)
        s_z = self.operators["s_z"]
        s_p = self.operators["s_p"]
        s_m = self.operators["s_m"]
        # set the matrix elements different from zero to the right values
        s_x[0, 1] = 0.5
        s_x[1, 0] = 0.5
        s_y[0, 1] = 1j*(-0.5)
        s_y[1, 0] = 1j*0.5
        s_z[0, 0] = -0.5
        s_z[1, 1] = 0.5
        s_p[1, 0] = 1.0
        s_m[0, 1] = 1.0
        # add the states
        self.add_state("spin_up")
        self.add_state("spin_down")
        self.add_state("empty")
        self.add_state("occupied")
        # for clarity
        state_up = self.states["spin_up"]
        state_down = self.states["spin_down"]
        state_empty = self.states["empty"]
        state_occupied = self.states["occupied"]
        # set the list elements different from zero to the right values
        state_up[1] = 1.0
        state_down[0] = 1.0
        state_occupied[1] = 1.0
        state_empty[0] = 1.0



class ElectronicSite(Site):
    """A site for electronic models
    
    You use this site for models where the single sites are electron
    sites. The Hilbert space is ordered such as:

    - the first state, labelled 0,  is the empty site,
    - the second, labelled 1, is spin down, 
    - the third, labelled 2, is spin up, and 
    - the fourth, labelled 3, is double occupancy.
    
    Notes
    -----
    Postcond: The site has already built-in the spin operators for: 

    - c_up : destroys an spin up electron,
    - c_up_dag, creates an spin up electron,
    - c_down, destroys an spin down electron,
    - c_down_dag, creates an spin down electron,
    - s_z, component z of spin,
    - s_p, raises the component z of spin,
    - s_m, lowers the component z of spin,
    - n_up, number of electrons with spin up,
    - n_down, number of electrons with spin down,
    - n, number of electrons, i.e. n_up+n_down, and
    - u, number of double occupancies, i.e. n_up*n_down.

    """
    def __init__(self):
        super(ElectronicSite, self).__init__(4)
        # add the operators
        self.add_operator("c_up")
        self.add_operator("c_up_dag")
        self.add_operator("c_down")
        self.add_operator("c_down_dag")
        self.add_operator("s_z")
        self.add_operator("n_up")
        self.add_operator("n_down")
        self.add_operator("n")
        self.add_operator("u")
        # for clarity
        c_up = self.operators["c_up"]
        c_up_dag = self.operators["c_up_dag"]
        c_down = self.operators["c_down"]
        c_down_dag = self.operators["c_down_dag"]
        s_z = self.operators["s_z"]
        n_up = self.operators["n_up"]
        n_down = self.operators["n_down"]
        n = self.operators["n"]
        u = self.operators["u"]
        # set the matrix elements different from zero to the right values
        c_up[0,2] = 1.0
        c_up[1,3] = 1.0
        c_up_dag[2,0] = 1.0
        c_up_dag[3,1] = 1.0
        c_down[0,1] = 1.0
        c_down[2,3] = 1.0
        c_down_dag[1,0] = 1.0
        c_down_dag[3,2] = 1.0
        s_z[1,1] = -1.0
        s_z[2,2] = 1.0
        n_up[2,2] = 1.0
        n_up[3,3] = 1.0
        n_down[1,1] = 1.0
        n_down[3,3] = 1.0
        n[1,1] = 1.0
        n[2,2] = 1.0
        n[3,3] = 2.0
        u[3,3] = 1.0
        # add the states
        self.add_state("empty")
        self.add_state("spin_down")
        self.add_state("spin_up")
        self.add_state("double")
        # for clarity
        state_empty = self.states["empty"]
        state_down = self.states["spin_down"]
        state_up = self.states["spin_up"]
        state_double = self.states["double"]
        # set the list elements different from zero to the right values
        state_empty[0] = 1.0
        state_down[1] = 1.0
        state_up[2] = 1.0
        state_double[3] = 1.0
