from moha.system.molecule import periodic
from moha.system.molecule.atom import Atom

__all__ = []

def load(filename):
    """Load a molecular geometry from a .xyz file.

    Parameters
    ----------
    filename: str
        The file to load the geometry from

    Returns
    -------
    atoms: list
        A list of atom.
    """
    if not isinstance(filename, str):
        raise TypeError("filename must be str")
    if not filename.endswith('.xyz'):
        raise ValueError("file format must be xyz")
    atoms = []
    with open(filename) as handle:
        size = int(next(handle))
        title = next(handle).strip()
        for _ in range(size):
            row = next(handle).split()
            tag = row[0]
            if isinstance(tag,int)  and 0<tag<121:
                element = list(periodic.table.values())[tag]
            elif tag.isdigit() and 0<int(tag)<121:
                element = list(periodic.table.values())[int(tag)]
            elif tag.lower() in periodic.table.keys():
                element = periodic.table[tag.lower()]
            else:
                raise ValueError("Element must be in the periodic table")
            coordinate = []
            for j in range(3):
                coordinate.append(float(row[j+1]))
            atom = Atom(element,coordinate)
            atoms.append(atom)
        handle.close()

    return atoms
