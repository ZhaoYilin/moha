import numpy as np
from moha.system.periodic import Element


class Atom(object):
    """
    Class for atom
    """
    def __init__(self,element,coordinate):
        for name, value in element.__dict__.items():
            setattr(self,name,value)
        self.coordinate = coordinate

class Molecule(object):
    """
    Molecular Class
        title
            type: string
            title of the system
        size
            type: intager
            number of atoms
        symmetry
            type: string
            point group of the molecule

    """
    def __init__(self,title,size,atoms=[],symmetry='C1'):
        self.title = title
        self.size = size
        self.atoms = atoms
        self.symmetry = symmetry

    def bond_length(self,i,j):
        """
        Calculate the the distances between atom i and j
        """
        A = np.array(self.atoms[i].coordinate)
        B = np.array(self.atoms[j].coordinate)
        return np.linalg.norm(A-B)

    def bond_angle(self,A,B,C):
        """
        Calculate the bond angle between atoms A B and C,
        where B is the central atom
        """
        return np.linalg.norm(A-B)

    def out_of_plane_angle(self,A,B,C,D):
        """
        Calculate out of plane angle, which is for atom A
        out of the plane consist by atoms B C D 
        """
        return np.linalg.norm(A-B)

    def dihedral_angle(self,A,B,C,D):
        """
        Calculate dihedral angle for atom connectivit A B 
        C D
        """
        return np.linalg.norm(A-B)

    def center_of_mass(self):
        """
        Calculate the center of mass for molecular
        """
        pass

