import numpy as np
from moha.log import log, timer
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




class MPNSolver(object):
    def __init__(self,energy):
        self.energy = energy

    @classmethod
    def mp2(cls,hfwavefunction,hamiltonian):
        occ = hfwavefunction.occ
        C = hfwavefunction.coefficient
        eorbitals = hfwavefunction.eorbitals
        Emp2 = 0.0
        if occ['alpha'] == occ['beta']:
            Eri = hamiltonian.operators['electron_repulsion'].basis_transformation(C)
            for i in range(occ['alpha']):
                for j in range(occ['alpha']):
                    for a in range(occ['alpha'],hamiltonian.dim):
                        for b in range(occ['alpha'],hamiltonian.dim):
                            Emp2 += Eri[i,a,j,b]*(2*Eri[i,a,j,b]-Eri[i,b,j,a])/(eorbitals[i] + eorbitals[j] -eorbitals[a] - eorbitals[b])

        elif occ['alpha'] != occ['beta']:
            for spin in C:
                Eri = hamiltonian.operators['electron_repulsion'].basis_transformation(C[spin])
                for i in range(occ[spin]):
                    for j in range(occ[spin]):
                        for a in range(occ[spin],hamiltonian.dim):
                            for b in range(occ[spin],hamiltonian.dim):
                                Emp2 += Eri[i,a,j,b]*(Eri[i,a,j,b]-0.5*Eri[i,b,j,a])/(eorbitals[spin][i] + eorbitals[spin][j] -eorbitals[spin][a] - eorbitals[spin][b])
        log.hline()
        log('{0:2s} {1:3f}'.format('Escf', hfwavefunction.Etot))
        log('{0:2s} {1:3f}'.format('Emp2', Emp2))
        log('{0:2s} {1:3f}'.format('Etot', hfwavefunction.Etot+Emp2))
        log.hline()

        return Emp2

    @classmethod
    def mp3(cls,hfwavefunction,hamiltonian):
        occ = hfwavefunction.occ
        C = hfwavefunction.coefficient
        pass
