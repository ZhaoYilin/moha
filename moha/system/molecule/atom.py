from moha.system.molecule.periodic import Element
from moha.system.basis import nwchem

__all__ = []
# Periodic table used in MoHa
# Default unit using in the peridoic table
# mass: amu
# radius: pm

class Atom(object):
    """Class for atom.

    Attributes
    ----------
    element : Element
        Element of the atom. 
        
    coordinate : 3-list
        Coordinate of the atom.

    Methods
    -------
    __init__(self,element,coordinate)
        Initialize the instance.
    
    assign_coordinate(self,coordinate)
        Assign coordinate to the instance.

    """
    def __init__(self, element: Element, coordinate: list):
        """Initialize the instance.

        Parameters
        ----------
        element : Element instance
            Chemical element instance.
    
        coordinate : 3-list
            Coordinate of the atom.
        """
        for name, value in element.__dict__.items():
            setattr(self,name,value)
        self.assign_coordinate(coordinate)
    
    def assign_coordinate(self,coordinate):
        """Assign coordinate to the instance.

        Parameters
        ----------
        coordinate : 3-list
            Coordinate of the atom.
        """
        if not isinstance(coordinate, list):
            raise TypeError("coordinate must be list") 
        if not len(coordinate)==3:
            raise ValueError("length of the coordinate must be 3") 
        self.coordinate = coordinate
