from moha.system.periodic import load_periodic
from moha.system.atom import Atom

import numpy as np

__all__ = ['Molecule']

class Molecule(object):
    """Molecular Class.

    Attributes
    ----------
    title : str
        Title of the system.

    size : int
        Number of atoms.

    symmetry : str
        Point group of the molecule.

    Property
    --------
    center_of_mass(self)
        Calculate the center of mass for molecular.

    Methods
    -------
    bond_length(self,i,j)
        
    """
    def __init__(self,title,size,atoms=[],symmetry='C1'):
        """Initialize the instance.

        Parameters
        ----------
        title : str
            Title of the system.

        size : int
            Number of atoms.

        symmetry : str
            Point group of the molecule.
        """
        self.title = title
        self.size = size
        self.atoms = atoms
        self.symmetry = symmetry

    @property
    def center_of_mass(self):
        """Calculate the center of mass for molecular.
        """
        pass
    
    @classmethod
    def from_file(cls,filename):
        """Generate N electron basis set.

        Parameters
        ----------
        truncation : list
            A list of truncation degree.

        Raises
        ------
        TypeError 
            If truncation is not a list.
        """
        periodic = load_periodic()
        #read molecule
        with open(filename) as f:
            size = int(next(f))
            title = next(f).strip()
            molecule = cls(title,size)
            for _ in range(size):
                row = next(f).split()
                tag = row[0]
                element = periodic[tag]
                coordinate = []
                for j in range(3):
                    coordinate.append(float(row[j+1]))
                atom = Atom(element,coordinate)

                molecule.atoms.append(atom)
            f.close()

        return molecule
        
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
            if item >= self.size:
                raise ValueError("Index of atom must be smaller than number of atoms")
        A = np.array(self.atoms[i].coordinate)
        B = np.array(self.atoms[j].coordinate)
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
            if item >= self.size:
                raise ValueError("Index of atom must be smaller than number of atoms")
        A = np.array(self.atoms[i].coordinate)
        B = np.array(self.atoms[j].coordinate)
        C = np.array(self.atoms[k].coordinate)
        return np.linalg.norm(A-B)

    def out_of_plane_angle(self,i,j,k,l):
        """Calculate out of plane angle, which is for atom A
           out of the plane consist by atoms B C D 

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
            if item >= self.size:
                raise ValueError("Index of atom must be smaller than number of atoms")
        A = np.array(self.atoms[i].coordinate)
        B = np.array(self.atoms[j].coordinate)
        C = np.array(self.atoms[k].coordinate)
        D = np.array(self.atoms[l].coordinate)

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
            if item >= self.size:
                raise ValueError("Index of atom must be smaller than number of atoms")
        A = np.array(self.atoms[i].coordinate)
        B = np.array(self.atoms[j].coordinate)
        C = np.array(self.atoms[k].coordinate)
        D = np.array(self.atoms[l].coordinate)

