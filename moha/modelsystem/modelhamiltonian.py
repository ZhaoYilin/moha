import numpy as np
from moha.modelsystem.auxiliary import tensor_enumerate
from moha.modelsystem.auxiliary import is_interaction
from moha.modelsystem.auxiliary import expectation

class ModelHamiltonian(object):
    """Base class for the Physical Model Hamiltonians: 
    1-D Hubbard, PPP, Ising, etc.
    """

    def __init__(self,lattice,site,pbc=True):
        """Initialize a ModelHamiltonian instance.

        Parameters
        ----------
        lattice : Lattice
            The lattice of the model hamiltonian.

        site : Site
            The site of the model hamiltonian

        nbasis : int
            The number of sites.

        pdb : bool
            Periodic boundary conditions. Default, pdb=true
        """
        self.lattice = lattice
        self.site = site
        self._pbc = pbc
        self._nbasis = lattice.Nsites
        self.operators = []
        self.one_body_matrix = np.zeros((self.nbasis,self.nbasis))
        self.two_body_tensor = np.zeros((self.nbasis,self.nbasis,self.nbasis,self.nbasis))

    @property
    def pbc(self):
        """The periodic boundary condition."""
        return self._pbc

    @property
    def nbasis(self):
        """The number of sites."""
        return self._nbasis

    def add_operator(self, operator):
        """Add operator to the model Hamiltonian."""
        self.operators.append(operator)
        self.compute_operator(operator)

    def compute_operator(self, operator):
        """Calculate the one-body term of the model Hamiltonian."""
        for i,listi in tensor_enumerate(self.lattice.sites):
            for j,listj in tensor_enumerate(self.lattice.sites):
                if is_interaction(self.lattice.shape,listi,listj,operator,False):
                    if operator.__class__.__name__ == 'OneBodyTerm':
                        if i != j:
                            self.one_body_matrix[i,j] = expectation(operator,self.site)
                            self.one_body_matrix[j,i] = expectation(operator,self.site)
                        elif i == j:
                            self.one_body_matrix[i,j] = expectation(operator,self.site)
                    elif operator.__class__.__name__ == 'TwoBodyTerm':
                        if i != j:
                            self.two_body_tensor[i,i,j,j] = expectation(operator,self.site)
                            self.two_body_tensor[j,j,i,i] = expectation(operator,self.site)
                        elif i == j:
                            self.two_body_tensor[i,i,j,j] = expectation(operator,self.site)

        if self._pbc == True:                
            for i,listi in tensor_enumerate(self.lattice.sites):
                for j,listj in tensor_enumerate(self.lattice.sites):
                    if is_interaction(self.lattice.shape,listi,listj,operator,True):
                        if operator.__class__.__name__ == 'OneBodyTerm':
                            if i != j:
                                self.one_body_matrix[i,j] = expectation(operator,self.site)
                                self.one_body_matrix[j,i] = expectation(operator,self.site)
                            elif i == j:
                                self.one_body_matrix[i,j] = expectation(operator,self.site)
                        elif operator.__class__.__name__ == 'TwoBodyTerm':
                            if i != j:
                                self.two_body_tensor[i,i,j,j] = expectation(operator,self.site)
                                self.two_body_tensor[j,j,i,i] = expectation(operator,self.site)
                            elif i == j:
                                self.two_body_tensor[i,i,j,j] = expectation(operator,self.site)

    def compute_overlap(self):
        """Calculate overlap of the model Hamiltonian, (identity matrix)."""
        return np.identity(self.nbasis)

    def build(self):
        ham = Hamiltonian(self._nbasis)
        overlap = self.compute_overlap()
        nuclear_attraction = np.zeros((self._nbasis,self._nbasis))
        Eri_list = []
        for i in range(self.nbasis):
            for j in range(self.nbasis)[0:i+1]:
                for k in range(self.nbasis)[0:i+1]:
                    for l in range(self.nbasis)[0:k+1]:
                        ij = i*(i+1)*0.5+j
                        kl = k*(k+1)*0.5+l
                        if ij>=kl:
                            Eri_list.append(self.two_body_tensor[i,j,k,l])

        ham.add_operator_by_value('nuclear_repulsion',1,0.)
        ham.add_operator_by_value('overlap',1,overlap)
        ham.add_operator_by_value('kinetic',1,self.one_body_matrix)
        ham.add_operator_by_value('nuclear_attraction',1,nuclear_attraction)
        ham.add_operator_by_value('electron_repulsion',2,Eri_list)
        return ham
