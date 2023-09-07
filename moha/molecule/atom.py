from moha.molecule.periodic import Element
from moha.basis import nwchem

__all__ = ['Atom']
# Periodic table used in MoHa
# Default unit using in the peridoic table
# mass: amu
# radius: pm


class Atom(object):
    
    def __init__(self, element: Element, coordinate: list):
        """_summary_

        Parameters
        ----------
        element : Element
            _description_
        coordinate : list
            _description_
        """
        for name, value in element.__dict__.items():
            setattr(self, name, value)
        self.assign_coordinate(coordinate)

    def __eq__(self, other, error=1e-5):
        """ 
        """
        if self.number != other.number:
            return False
        if np.allclose(self.coordinate,other.coordinate,error):
            return False
        else:
            return True

    def assign_coordinate(self, coordinate: list):
        """_summary_

        Parameters
        ----------
        coordinate : list
            _description_

        Raises
        ------
        TypeError
            _description_
        ValueError
            _description_
        """
        if not isinstance(coordinate, list):
            raise TypeError("coordinate must be list")
        if not len(coordinate) == 3:
            raise ValueError("length of the coordinate must be 3")
        self.coordinate = coordinate
