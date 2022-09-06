from moha.system.molecule import periodic
from moha.system.molecule.atom import Atom

__all__ = []

def load(geometry):
    """Load molecular geometry from a python iterator.

    Parameters
    ----------
    geometry: list,tuple
        The geometry of molecular.

    Returns
    -------
    atoms: list
        A list of atom.
    """
    if not isinstance(geometry, (list,tuple)):
        raise TypeError("geometry must be iterator")
    atoms = []

    for row in geometry:
        tag = row[0]
        if isinstance(tag,int)  and 0<tag<121:
            element = list(periodic.table.values())[tag]
        elif tag.isdigit() and 0<int(tag)<121:
            element = list(periodic.table.values())[int(tag)]
        elif tag.lower() in periodic.table.keys():
            element = periodic.table[tag.lower()]
        else:
            raise ValueError("Element must be in the periodic table")
        coordinate = list(row[1:])
        atom = Atom(element,coordinate)
        atoms.append(atom)

    return atoms
