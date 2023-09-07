from moha.system.wavefunction.ci_wavefunction import CIWaveFunction
from moha.system.basis.basis import SlaterDeterminant
from moha.system.basis.basis_set import FockSpace
from moha.io.log import log,timer

import numpy as np

__all__ = []

def cisd(ham,wfn):
    """Kernel of the solver.

    Returns
    -------
    results : dict
        CISD calculation results.
    """
    log.hline(char='=')
    log.blank()
    log('CISD Section'.format())

    nspatial = 7
    nelec = 10
    nspin = nspatial*2

    reference = wfn.ansatz
    cibs = FockSpace.ci_space(reference,2)
    
    C = wfn.coefficient

    Hcore = ham['Hcore'].basis_transform(C,C)
    Hcore = np.repeat(Hcore,2,axis=0)
    h1e = np.repeat(Hcore,2,axis=1)
    # Make H block diagonal
    spin_ind = np.arange(h1e.shape[0], dtype=np.int) % 2
    h1e *= (spin_ind.reshape(-1, 1) == spin_ind)

    Eri = ham['Eri'].basis_transform(C,C,C,C)

    g2e = Eri.double_bar
    ham.update({'h1e':h1e})
    ham.update({'g2e':g2e})
    print(len(cibs))
    ci_matrix = np.zeros((len(cibs),len(cibs)))
    for i in range(len(cibs)):
        for j in range(i+1):
            ci_matrix[i,j] = ham(cibs[i],cibs[j])
            ci_matrix[j,i] = ci_matrix[i,j]

    e_ci, vec_ci = np.linalg.eigh(ci_matrix)
    Eci = e_ci[0]
    
    #Build the ci wavefunction.
    #ci_wfn = CIWaveFunction(nelec,nspatial,{},cibs,vec_ci[0])
    
    log.hline()
    log('CISD Results')
    log.hline()
    log('{0:<20s}{1:<20.6f}'.format('Electronic energy',Eci))
    log('{0:<20s}{1:<20.6f}'.format('Electronic energy',Eci+float(ham['Nri'])))
    log.blank()
    
    results = {
    "success":True,
    "Eci":Eci,
    }

    return results
