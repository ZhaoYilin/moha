import os
from moha.system.periodic import load_periodic
from moha.system.atom import Atom
from moha.system.molecule import Molecule

def load_xyz(filename):
    """Load a molecular geometry from a .xyz file.

       **Argument:**

       filename
            The file to load the geometry from

       **Returns:** dictionary with ``title`, ``coordinates`` and ``numbers``.
    """
    periodic = load_periodic()
    #read molecule
    with open(filename) as f:
        size = int(next(f))
        title = next(f).strip()
        molecule = Molecule(title,size)
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


