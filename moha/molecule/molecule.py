import itertools
import math
import numpy as np

from moha.molecule import periodic
from moha.molecule.atom import Atom
from moha.basis.basis_set import BasisSet
from moha.io.log import log

__all__ = ['Molecule']

class Molecule(list):
    """Molecule Class.

    Attributes
    ----------
    pg: bool or str
        Point group symmetry of the molecule.

    Property
    --------
    E_nuc(self)
        Nucelar repulsion energy.

    center_of_mass(self)
        Center of mass for molecular.

    Methods
    -------
    assign_point_group(self,pg)
        Assign point group symmetry to the molecule.

    bond_length(self,i,j)
        Calculate the the distances between atom i and j.

    bond_angle(self,i,j,k)
        Calculate the bond angle between atoms i,j,k.
 
    dihedral_angle(self,i,j,k,l)
        Calculate dihedral angle for atom connectivit i,j,k,l.
    """
    def __new__(cls):
        """Generate new molecule object.
        """
        obj = list().__new__(cls)

        return obj

    def __eq__(self,other,error=1e-3):
        """
        """
        geo_A = self.geometry
        geo_B = other.geometry
        geo_A = np.sort(geo_A,axis=0)
        geo_B = np.sort(geo_B,axis=0)
        return np.allclose(geo_A,geo_B,atol=error)

    @classmethod
    def build(cls,geometry,pg=False):
        """Generate N electron basis set.

        Parameters
        ----------
        geometry : str or iterator
            Geometry object.

        pg: bool or str
            Point group symmetry of the molecule.

        Raises
        ------
        TypeError 
            If geometry is unknow format.
        """
        log('Molecule Section'.format())
        log.hline()

        molecule = Molecule()
        molecule.assign_point_group(pg)
        if isinstance(geometry, str):
            if geometry.endswith('.xyz'):
                from moha.molecule import xyz
                atoms = xyz.load(geometry)

        elif isinstance(geometry,(list,tuple,dict)):
            from moha.molecule import iterator
            atoms = iterator.load(geometry)

        else:
            raise ValueError("Unknown format")

        log('{0:<20s} {1:<20s} {2:<20s} {3:<20s}'.format('Element', 'X', 'Y','Z'))
        log.hline()
        for atom in atoms:
            molecule.append(atom)
            log('{0:<20s} {1:<20f} {2:<20f} {3:<20f}'.format(atom.symbol, 
            atom.coordinate[0], atom.coordinate[1],atom.coordinate[2]))
        log.blank()


        return molecule

    def append(self,atom):
        """Add atom to the molecule.

        Parameters
        ----------
        atom : Atom
            Atom instance.
        
        Raises
        ------
        TypeError
            If atom is not an Atom instance.
        """
        if not isinstance(atom, Atom):
            raise TypeError("Atom must be a Atom instance")
        super().append(atom)

    def assign_point_group(self,pg):
        """Assign point group symmetry to the molecule.

        Parameters
        ----------
        pg: bool or str
            Point group symmetry of the molecule.

        Raises
        ------
        TypeError
            If pg is not bool or str.
        """
        if not isinstance(pg, bool or str):
            raise TypeError("Parameters pg must be bool or str.")
        self.pg = pg
                                
    @property
    def e_nuc(self):
        """Nucelar repulsion energy.

        Return
        ------
        energy: float
            Nuclear repulsion energy.
        """
        energy = 0.
        for i,A in enumerate(self):
            for j,B in enumerate(self):
                if j>i:
                    energy += (A.number*B.number)/self.bond_length(i,j)
        return energy

    def bond_length(self,i,j):
        """Calculate the the distances between atom i and j.

        Parameters
        ----------
        i : int
            Index of first atom.
        
        j : int
            Index of second atom.
            
        Raises
        ------
        TypeError
            If index of atom is not int.
        ValueError
            If index of atom is not zero or positive number.            
            If index of atom is not samller than number of atoms.            
        """
        for item in (i,j):
            if not isinstance(item,int):
                raise TypeError("Index of atom must be int")
            if item < 0:
                raise ValueError("Index of atom must be none negative number")
            if item >= len(self):
                raise ValueError("Index of atom must be smaller than number of atoms")
        A = np.array(self[i].coordinate)
        B = np.array(self[j].coordinate)
        return np.linalg.norm(A-B)

    def bond_angle(self,i,j,k):
        """Calculate the bond angle between atoms i,j,k.
           Where atom j is the central atom.

        Parameters
        ----------
        i : int
            Index of first atom.
        
        j : int
            Index of second atom.
            
        k : int
            Index of third atom.

        Raises
        ------
        TypeError
            If index of atom is not int.
        ValueError
            If index of atom is not zero or positive number.            
            If index of atom is not samller than number of atoms.            
        """
        for item in (i,j,k):
            if not isinstance(item,int):
                raise TypeError("Index of atom must be int")
            if item < 0:
                raise ValueError("Index of atom must be none negative number")
            if item >= len(self):
                raise ValueError("Index of atom must be smaller than number of atoms")
        A = np.array(self[i].coordinate)
        B = np.array(self[j].coordinate)
        C = np.array(self[k].coordinate)
        AB = A-B
        BC = B-C
        
        inner = np.inner(AB, BC)
        norms = np.linalg.norm(AB) * np.linalg.norm(BC)
        cos = inner / norms
        rad = np.arccos(np.clip(cos, -1.0, 1.0))
        return rad

    def dihedral_angle(self,A,B,C,D):
        """Calculate dihedral angle for atom connectivit A,B,C,D.

        Parameters
        ----------
        i : int
            Index of first atom.
        
        j : int
            Index of second atom.
            
        k : int
            Index of third atom.

        l : int
            Index of fourth atom.

        Raises
        ------
        TypeError
            If index of atom is not int.
        ValueError
            If index of atom is not zero or positive number.            
            If index of atom is not samller than number of atoms.            
        """
        for item in (i,j,k):
            if not isinstance(item,int):
                raise TypeError("Index of atom must be int")
            if item < 0:
                raise ValueError("Index of atom must be none negative number")
            if item >= len(self):
                raise ValueError("Index of atom must be smaller than number of atoms")
        A = np.array(self[i].coordinate)
        B = np.array(self[j].coordinate)
        C = np.array(self[k].coordinate)
        D = np.array(self[l].coordinate)

    def normal_vector(self,i,j):
        """
        """
        r_i = np.array(self[i].coordinate)
        r_j = np.array(self[j].coordinate)
        n = (r_i-r_j)/np.linalg.norm(r_i-r_j)
        return n

    @property
    def geometry(self):
        """  This function calculates the center of mass of a given data array.
        The format of the data array is a two dimensional array, which contains
        the x, y, z Cartesion coordinates and the corresponding mass values for
        each x, y, z coordinates.

        It returns a 1-D data array.
        """
        geometry = []
        for atom in self:
            single_line = [atom.number] + atom.coordinate
            geometry.append(single_line)
        return geometry 

    @property
    def center_of_mass(self):
        """  This function calculates the center of mass of a given data array.
        The format of the data array is a two dimensional array, which contains
        the x, y, z Cartesion coordinates and the corresponding mass values for
        each x, y, z coordinates.

        It returns a 1-D data array.
        """
        geometry = []
        mass = [] 
        for atom in self:
            geometry.append(atom.coordinate)
            mass.append(atom.mass)
        return np.average(geometry, axis=0, weights=mass)

    @property
    def moment_of_inertia(self):
        """  This function calculates the Elements of inertia tensor for the
        center-moved coordinates.
    
        The structure of the array is a two dimensional array, which contains
        the mass of elements (the first column) and their corresponded x, y, z
        Cartesion coordinates.

        Return
        ------
        I : ndarray((3,3))
            The inertia matrix.
        """
        geometry = []
        mass = [] 
        for atom in self:
            geometry.append(atom.coordinate)
            mass.append(atom.mass)
        geometry = np.array(geometry)
        mass = np.array(mass)
        com_coord = np.average(geometry, axis=0, weights=mass)
        new_coord = geometry - com_coord
        I_xx = (mass * np.sum(np.square(new_coord[:,1:3:1]),axis=1)).sum()
        I_yy = (mass * np.sum(np.square(new_coord[:,0:3:2]),axis=1)).sum()
        I_zz = (mass * np.sum(np.square(new_coord[:,0:2:1]),axis=1)).sum()
        I_xy = (-1 * mass * np.prod(new_coord[:,0:2:1],axis=1)).sum()
        I_yz = (-1 * mass * np.prod(new_coord[:,1:3:1],axis=1)).sum()
        I_xz = (-1 * mass * np.prod(new_coord[:,0:3:2],axis=1)).sum()
        I = np.array([[I_xx, I_xy, I_xz],
		    [I_xy, I_yy, I_yz],
		    [I_xz, I_yz, I_zz]])
        return I

    @property
    def principal_moments_of_inertia(self):
        """
        """
        I = self.moment_of_inertia
        w, v = np.linalg.eigh(I)
        return w, v

    @property
    def distance_matrix(self):
        """interatomic distance matrix
        """
        natom = len(self)
        D = np.zeros((natom,natom))
        for i, j in itertools.product(range(natom), repeat=2):
            if i <= j:
                D[i,j] = self.bond_length(i,j)
        return D + D.T - np.diag(np.diag(D))

    @property
    def equivalent_atoms(self):
        """
        """
        EA = []
        element_set = set([atom.number for atom in self])
        for element in element_set:
            sub_mol = self.__class__()   
            for atom in self:
                if atom.number==element:
                    sub_mol.append(atom)
            EA.append(sub_mol)
        return EA
        
    def symmetrically_equivalent_atoms(self,error=1e-3):
        """SEA
        """
        error_order = -int(math.log10(error))

        SEA = []
        EA = self.equivalent_atoms
        for sub_mol in EA:
            if len(sub_mol)==1:
                SEA.append(sub_mol)
            else:
                D_sub = sub_mol.distance_matrix
                D_sub = np.around(D_sub,error_order)
                D_sub = np.sort(D_sub)
                distance_set = np.unique(D_sub, axis=0)
                for i,distance_i in enumerate(distance_set):
                    subsub_mol = self.__class__()   
                    for j,distance_j in enumerate(D_sub):
                        if np.allclose(distance_i,distance_j):
                            subsub_mol.append(sub_mol[j])
                SEA.append(subsub_mol)
        return SEA

    def translation(self, displacement):
        """ Translate the molecule by the given displacement.

        Parameters
        ----------
        displacement : List[float, float, float]
            Vector of the displacement.

        Return
        ------
        new_mol : Molecule
            New molecule.
        """
        geometry = np.array(self.geometry)
        geometry[:,1:] -= displacement
        new_mol = self.build(geometry.tolist(),self.pg)
        return new_mol

    def basis_transformation(self, basis_set):
        """ Set molecular orientation along new basis set.

        Parameters
        ----------
        basis_set : ndarray((3,3))
            Stack of three unitary orthogonal vectors.

        Return
        ------
        new_mol : Molecule
            New molecule.
        """
        geometry = np.array(self.geometry)
        geometry[:,1:] = np.dot(geometry[:,1:],basis_set.T)
        new_mol = self.build(geometry.tolist(),self.pg)
        return new_mol

    def standary_orientation(self, error=1e-3):
        """ Set molecular orientation to x,y,z axes. Here z axis is the principal axis.

        Return
        ------
        new_mol : Molecule
            New molecule with standary orientation.
        """
        error_order = -int(math.log10(error))
        v, w = self.principal_moments_of_inertia
        v = np.around(v,error_order)
        w = w.T
        # Move principal axis to z axis.
        if v[0]==0 and v[1]==v[2]:
            # Linear
            standary_basis = np.array([w[1],w[2],w[0]])
        elif v[0]==v[1]==v[2]:
            # Spherical Top
            standary_basis = np.eye(3)
        elif v[0]==v[1]<v[2]:
            # Oblate symmetric Top
            standary_basis = w
        elif v[0]<v[1]==v[2]:
            # Prolate symmetric Top
            standary_basis = np.array([w[1],w[2],w[0]])
        elif v[0]!=v[1]!=v[2]:
            # Asymmetric Top
            standary_basis = np.array([w[2],w[0],w[1]])
        else:
            return False

        new_mol = self.basis_transformation(standary_basis)
        return new_mol

