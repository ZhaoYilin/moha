import os
import numpy as np
from moha.system.basis import GaussianOrbital,BasisSet
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


def load_nwchem(symbols,coordinates,filename):
    """
    read basis set
    """
    basis_set = BasisSet()
    #read the file and store them as list
    path = os.path.dirname(__file__) 
    with open(path + "/basis/"+filename) as f:
        content=[]
        for line in f:
            line = line[:line.find('#')].strip()
            if len(line) == 0 or line.startswith('BASIS'):
                continue
            content.append(line)
        f.close()

    for i,symbol in enumerate(symbols):
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
        for j,location in enumerate(locations[:-1]):
            words = content[location].split()
            for k,shell_symbol in enumerate(list(words[1].lower())):
                shells = shell_expansion(shell_symbol)
                for shell in shells:
                    exps = []
                    coefs = []
                    for l in range(locations[j]+1,locations[j+1]):
                        words = content[l].split()
                        exps.append(float(words[0]))
                        coefs.append(float(words[k+1]))
                    orb = GaussianOrbital.spatial(i,coordinates[i],shell,exps,coefs)
                    basis_set.size +=1
                    basis_set.basis.append(orb)

    return basis_set
