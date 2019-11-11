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



class ExcitationEnergy(object):
    
    @classmethod
    def tdhf(cls,hfwavefunction,hamiltonian):
        occ = hfwavefunction.occ
        Nelec = occ['alpha'] + occ['beta']
        dim = hamiltonian.dim
        C = hfwavefunction.coefficient
        Nov = Nelec*(2*dim - Nelec)
        #Transfer Fock integral from spatial to spin basis
        fs = spinfock(hfwavefunction.eorbitals)
        #Transfer electron repulsion integral from atomic basis
        #to molecular basis
        hamiltonian.operators['electron_repulsion'].basis_transformation(C)
        #build double bar integral <ij||kl>
        spinints = hamiltonian.operators['electron_repulsion'].double_bar
        
        A = np.zeros((Nov,Nov))
        B = np.zeros((Nov,Nov))
        I = -1
        for i in range(0,Nelec):
            for a in range(Nelec,dim*2):
                I += 1
                J = -1
                for j in range(0,Nelec):
                    for b in range(Nelec,dim*2):
                        J += 1
                        A[I,J] = (fs[a,a] - fs[i,i])*(i==j)*(a==b)+spinints[a,j,i,b]
                        B[I,J] = spinints[a,b,i,j]
        M = np.bmat([[A,B],[-B,-A]])
        ETD,CTD = np.linalg.eig(M)
        print(ETD)

    @classmethod
    def cis(cls,hfwavefunction,hamiltonian):
        occ = hfwavefunction.occ
        Nelec = occ['alpha'] + occ['beta']
        dim = hamiltonian.dim
        C = hfwavefunction.coefficient
        Nov = Nelec*(2*dim - Nelec)
        #Transfer Fock integral from spatial to spin basis
        fs = spinfock(hfwavefunction.eorbitals)
        #Transfer electron repulsion integral from atomic basis
        #to molecular basis
        hamiltonian.operators['electron_repulsion'].basis_transformation(C)
        #build double bar integral <ij||kl>
        spinints = hamiltonian.operators['electron_repulsion'].double_bar
        
        H = np.zeros((Nov,Nov))
        I = -1
        for i in range(0,Nelec):
            for a in range(Nelec,dim*2):
                I += 1
                J = -1
                for j in range(0,Nelec):
                    for b in range(Nelec,dim*2):
                        J += 1
                        H[I,J] = (fs[a,a] - fs[i,i])*(i==j)*(a==b)+spinints[a,j,i,b]
        ECIS,CCIS = np.linalg.eig(H)
        print(ECIS)



