from moha.system.basis.basis import ContractedGaussian,Orbital
from moha.system.basis.basis_set import BasisSet
from moha.io.log import log

import os
import numpy as np

__all__ = []

def isDigit(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

def shell_expansion(shell):
    s_shell = [(0,0,0)]
    p_shell = [(1,0,0),(0,1,0),(0,0,1)]
    d_shell = [(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]
    shell_directory = {'s': s_shell, 'p': p_shell, 'd': d_shell}
    shells = shell_directory[shell]
    return shells


def load(atom,basis_set,filename):
    """Load the basis set.

    Parameters
    ----------
    coordinate: 3-list or 3-tuple
        Coordinate of the basis.

    filename: str
        Name of the basis set file.

    Return
    ------
    basis set: BasisSet
        Basis set instance.
    """
    if not isinstance(filename, str):
        raise TypeError("filename must be str.")
    if not filename.endswith('.nwchem'):
        raise ValueError("file format must be nwchem.")

    symbol = atom.symbol
    coordinate = atom.coordinate

    #read the file and store them as list
    path = os.path.dirname(__file__) 
    with open(path + "/database/"+filename) as handle:
        content=[]
        for line in handle:
            line = line[:line.find('#')].strip()
            if len(line) == 0 or line.startswith('BASIS'):
                continue
            content.append(line)
        handle.close()

    #read the location of all the line initial with symbol
    locations = []
    for index,line in enumerate(content):
        words = line.split()
        if words[0].lower()==symbol.lower():
            locations.append(index)
    for index,line in enumerate(content[locations[-1]+1:]):
        words = line.split()
        if isDigit(words[0]):
            continue
        else:
            locations.append(locations[-1]+index+1)
            break

    #build orbital and add it to basis set 
    for i,location in enumerate(locations[:-1]):
        words = content[location].split()
        principal_number = i
        for j,shell_symbol in enumerate(list(words[1].lower())):
            shells = shell_expansion(shell_symbol)
            for shell in shells:
                exps = []
                coefs = []
                for k in range(locations[i]+1,locations[i+1]):
                    words = content[k].split()
                    exps.append(float(words[0]))
                    coefs.append(float(words[j+1]))
                orb = Orbital()
                cg = ContractedGaussian(coordinate,principal_number,shell,coefs,exps)
                orb.append(cg)
                basis_set.append(orb)

    return basis_set
