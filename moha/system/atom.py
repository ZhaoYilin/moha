from moha.system.periodic import Element

__all__ = ['Atom']

class Atom(object):
    """Class for atom.

    Attributes
    ----------
    name : str
        Name of the atom.

    value : int
        Atomic number of the atom       
        
    coordinate : 3-list
        Coordinate of the atom.

    Methods
    -------
    __init__(self,element,coordinate)
        Initialize the instance.
    """
    def __init__(self,element,coordinate):
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
        self.coordinate = coordinate

