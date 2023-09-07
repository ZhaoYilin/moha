import numpy as np
import copy

from moha.hamiltonian.operator import *
from moha.basis.n_electrons import SlaterDeterminant
from moha.io.log import log,timer

__all__ = ['Hamiltonian']

class Hamiltonian(dict):
    """Chemical Hamiltonian for a Schrodinger equation.

    Methods
    -------
    __init__(self, nspatial)
        Initialize the Hamiltonian.

    assign_operators(self, operator)
        Assigns operator integral to the Hamiltonian.

    ClassMethods
    ------------
    build(cls,molecule,basis_set)
        Build the Chemical Hamiltonian
    """
    def __new__(cls):
        """Generate new molecule object.
        """
        obj = dict().__new__(cls)

        return obj

    def __call__(self,sd1,sd2):
        """Calculate matrix element between two determinants.

        Parameters
        ----------
        sd1 : SlaterDeterminant
            The first SlaterDeterminant instance.
        
        sd2 : SlaterDeterminant
            The second SlaterDeterminant instance.
        
        Raises
        ------
        TypeError
            If sd1 is not SlaterDeterminant instance.
            If sd2 is not SlaterDeterminant instance.

        Returns
        -------
        value : float
            The integral value of CI Hamiltonian.
        """
        if not isinstance(sd1,SlaterDeterminant):
            raise TypeError("Parameters det1 must be SlaterDeterminant instance.")
        if not isinstance(sd2,SlaterDeterminant):
            raise TypeError("Parameters det2 must be SlaterDeterminant instance.")
        degree = sd1.degree_of_excitation(sd2)

        #Calculate matrix element whose configurations are identical.
        if degree == 0:
            O1 = 0.0
            O2 = 0.0
            for i in sd1.occ_indices:
                O1 += self['h1e'][i,i]
            for i in sd1.occ_indices:
                for j in sd1.occ_indices:
                    O2 += self['g2e'][i,j,i,j]
            value = O1+0.5*O2
        #Calculate matrix element whose configuration differ by 1 spin orbitals.
        elif degree == 1:
            O1 = 0.0
            O2 = 0.0
            sign = sd1.phase(sd2)
            occ_indices = copy.deepcopy(sd1.occ_indices)
            a_list = sd1.annihilation_list(sd2)
            c_list = sd1.creation_list(sd2)

            for i in a_list:
                if i in occ_indices:
                    occ_indices.remove(i)
            for j in c_list:
                if j in occ_indices:
                    occ_indices.remove(j)
            O1 += self['h1e'][a_list[0],c_list[0]]
            for k in occ_indices:
                O2 += self['g2e'][a_list[0],k,c_list[0],k]
            value = O1+O2
            value *=sign
        #Calculate matrix element whose configuration differ by 2 spin orbitals.
        elif degree == 2:
            value = 0.0
            sign = sd1.phase(sd2)
            a_list = sd1.annihilation_list(sd2)
            c_list = sd1.creation_list(sd2)
            value += self['g2e'][a_list[0],a_list[1],c_list[0],c_list[1]]
            value *=sign
        else:
            value = 0.0

        return value

    @classmethod
    @timer.with_section('Integral')
    def build(cls,mol,orb):
        """Build the chemical hamiltonian.

        Parmeters
        ---------
        mol : Molecule
            Instance of Molecule.

        orb : BasisSet
            Instance of BasisSet.
        """
        #Section log: header.    
        log.hline(char='=')
        log.blank()
        log("Hamiltonian Section".format())
        log.hline()

        #Create new object instance.
        ham = cls()
        
        #Build operators.
        Nri = NuclearRepulsion.build(mol)
        S = Overlap.build(orb)
        T = Kinetic.build(orb)
        V = NuclearAttraction.build(mol,orb)
        Hcore = T+V
        Eri = ElectronRepulsion.build(orb)
        
        #Assign operators to Hamiltonian.
        ham['Nri'] = Nri
        ham['S'] = S
        ham['T'] = T
        ham['V'] = V
        ham['Hcore'] = Hcore
        ham['Eri'] = Eri
        
        #Section log: results.    
        log('{0:<20s} {1:<20s}'.format('Operator','Use'))
        log.hline()
        for key,value in ham.items():
            log('{0:<20s} ham[\'{1:s}\']'.format(key,key))
        log.blank()

        return ham
    
    @timer.with_section('Integral')
    def load(self,*tensors):
        """Build the chemical hamiltonian

        Parmeters
        ---------
        molecule : Molecule instance
            Instance of molecule class

        basis_set : Basis Set instance
            Instance of one electron basis set
        """
        import os
        from moha.system.operator.operator import Operator

        #Section log: header.    
        log.hline(char='=')
        log.blank()
        log("Hamiltonian Section".format())
        log.hline()

        #Load and assign operators.
        for tensor in tensors:
            path,fullname = os.path.split(tensor)
            name,extension = os.path.splitext(fullname)
            key,value = name,np.load(tensor)
            if value.ndim==2:
                operator = Operator(value)
            elif value.ndim==4:
                operator = ElectronRepulsion(value)
            else:
                raise ValueError("Dimension of operator out of range.")
            self.update({key:operator})

        #Section log: results.    
        log('{0:<20s} {1:<20s}'.format('Operator','Use'))
        log.hline()
        for key,value in self.items():
            log('{0:<20s} ham[\'{1:s}\']'.format(key,key))
        log.blank()

        return self

    def dump(self):
        """Build the chemical hamiltonian
        """
        for key,value in self.items():
            np.save(key,value)
