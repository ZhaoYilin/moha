import numpy as np

def spinfock(eorbitals):
    """
    """
    if type(eorbitals) is np.ndarray:
        dim = 2*len(eorbitals)
        fs = np.zeros(dim)
        for i in range(0,dim):
            fs[i] = eorbitals[i//2]
        fs = np.diag(fs) # put MO energies in diagonal array
    elif type(eorbitals) is dict:
        dim = 2*len(eorbitals['alpha'])
        fs = np.zeros(dim)
        for i in range(0,dim):
            if i%2==0:
                fs[i] = eorbitals['alpha'][i//2]
            elif i%2==0:
                fs[i] = eorbitals['beta'][i//2]
        fs = np.diag(fs) # put MO energies in diagonal array
    return fs
