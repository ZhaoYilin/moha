from moha.system.wavefunction.ci_wavefunction import CIWaveFunction
from moha.system.basis.basis import SlaterDeterminant
from moha.system.basis.basis_set import FockSpace
from moha.io.log import log,timer

import numpy as np

__all__ = []

def cis(ham,wfn):
    """Kernel of the solver.

    Returns
    -------
    results : dict
        CIS calculation results.
    """
    log.hline(char='=')
    log.blank()
    log('CIS Section'.format())

    nspatial = ham.nspatial
    nelec = wfn.symmetry.n
    nspin = nspatial*2

    reference = wfn.ansatz
    cibs = FockSpace.ci_space(reference,1)
    
    C = wfn.coefficient
    Hcore = ham['Hcore'].basis_transformation(C)
    h1e = Hcore.spin_orbital_basis_integral
    Eri = ham['Eri'].basis_transformation(C)
    g2e = Eri.double_bar
    ham.update({'h1e':h1e})
    ham.update({'g2e':g2e})

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
    log('CIS Results')
    log.hline()
    log('{0:<20s}{1:<20.6f}'.format('Electronic energy',Eci))
    log.blank()
    
    results = {
    "success":True,
    "Eci":Eci,
    }

    return results
