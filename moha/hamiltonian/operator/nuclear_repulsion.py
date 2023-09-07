from moha.hamiltonian.operator.base import Operator

class NuclearRepulsion(Operator):
    """Nuclear repulsion operator class.
    
    Methods
    -------
    """
    @classmethod
    def build(cls,mol):
        """Build the nuclear repulsion operator.

        Parameters
        ----------
        mol : Molecule 
            Molecule instance.

        Raises
        ------
        TypeError
            If parameter mol is not instance of Molecule.

        Return
        ------
        operator: NuclearRepulsion
            Return NuclearRepulsion operator instance.
        """
        e_nuc = mol.e_nuc
        operator = cls([e_nuc])

        return operator

    def basis_transform(self):
        """Rotate orbitals with a transformation matrix.

        Return
        ------
        NotImplemented
        """
        return NotImplemented
