import numpy as np

__all__ = []

class WaveFunction(object):
    """Base wavefunction class. Here list some of properties and symmetries of 
    the exact wave function as methods that may guide us in the construction of 
    models.

    Properties
    ----------
    n_electron(self)
        The exact state is a function of the space and spin coordinates of N
        electrons.    

    antisymmetric(self)
        The exact wave function is antisymmetric with respect to the 
        permutation of any pair of electrons.

    normalization(self)
        For a bound state, the eaxct wave function is square-integrable and 
        hence normalizable.

    variational(self)
        The exact wave function is variational in the sense that, for all 
        possible variations that are orthogonal to the wave function, the energy
        is stable.

    size_extensive(self)
        The exact wave function is size-extensive.

    spin(self)
        In nonrelativistic theory, the exact stationary states are 
        eigenfunctions of the total and projected spin operators.

    spatial_symmetry(self)
        Within the Born-Oppenheimer approximation, the exact stationary 
        states form a basis for an irreducible representation of the 
        molecular point group.

    electronic_coulomb_cusp_condition(self)
        Owing to the presence of the Coulomb optential, the molecular 
        electronic Hamiltonian becomes singular when two electrons coincide 
        in space.

    long_range_exponential_decay(self)
        At large distances, the electron density decays exponentially.

    gauge_transformation(self)
        The exact wave function changes in a characteristic way under gauge 
        transformations of the potentials associated with electromagnetic 
        fields.
    """
    @property
    def n_electron(self):
        """The exact state is a function of the space and spin coordinates of N
        electrons.    

        Return
        ------
        False
        """
        return False

    @property
    def antisymmetric(self):
        """The exact wave function is antisymmetric with respect to the 
        permutation of any pair of electrons.

        Return
        ------
        False
        """
        return False

    @property
    def normalization(self):
        """For a bound state, the eaxct wave function is square-integrable and 
        hence normalizable.

        Return
        ------
        False
        """
        return False
 
    @property
    def variational(self):
        """The exact wave function is variational in the sense that, for all 
        possible variations that are orthogonal to the wave function, the energy
        is stable.

        Return
        ------
        False
        """
        return False
        
    @property
    def size_extensive(self):
        """The exact wave function is size-extensive.

        Return
        ------
        False
        """
        return False
        
    @property
    def spin(self):
        """In nonrelativistic theory, the exact stationary states are 
        eigenfunctions of the total and projected spin operators.

        Return
        ------
        False
        """
        return False

    @property
    def spatial_symmetry(self):
        """Within the Born-Oppenheimer approximation, the exact stationary 
        states form a basis for an irreducible representation of the 
        molecular point group.

        Return
        ------
        False
        """
        return False
        
    @property
    def electronic_coulomb_cusp_condition(self):
        """Owing to the presence of the Coulomb optential, the molecular 
        electronic Hamiltonian becomes singular when two electrons coincide 
        in space.

        Return
        ------
        False
        """
        return False
 
    @property
    def long_range_exponential_decay(self):
        """At large distances, the electron density decays exponentially.

        Return
        ------
        False
        """
        return False
 
    @property
    def gauge_transformation(self):
        """The exact wave function changes in a characteristic way under gauge 
        transformations of the potentials associated with electromagnetic 
        fields.

        Return
        ------
        False
        """
        return False
