from moha.io.log import log, timer

import numpy as np

__all__ = ['MPSolver']

class MPSolver(object):
    """Moller-Plesset perturbation solver.
     
    Attributes
    ----------
    ham: Hamiltonian
        Chemical Hamiltonian.

    wfn: HFWaveFunction
        Hartree Fock wavefunction.
    
    Methods
    -------
    __init__(self,ham,wfn)
        Initialize the solver.

    mp2(self)
        MP2 solver.

    assign_hamiltonian(self,ham)
        Assign the chemical Hamiltonian to the solver.

    assign_wavefunction(self,wfn)
        Assign the Hartree Fock wavefunction to the solver.
    """
    def __init__(self,ham,wfn):
        """Initialize the solver.
 
        Parameters
        ----------
        ham: Hamiltonian
            Chemical Hamiltonian.

        wfn: HFWaveFunction
            Hartree Fock wavefunction.
        """
        self.assign_hamiltonian(ham)
        self.assign_wavefunction(wfn)        

    @timer.with_section('MP2')
    def mp2(self):
        """MP2 solver.

        Returns
        -------
        E: float
            MP2 energy.
        """
        from moha.mp.mp2 import mp2
        log.hline(char='=')
        log.blank()
        log('MP2 Section'.format())

        ham = self.ham
        wfn = self.wfn
        
        # Prepare parameters for mp2 algorithm
        Eri = ham['Eri']
        C = wfn.coefficient
        Emo = wfn.orbital_energies
        na = wfn.ansatz.nalphas
        nb = wfn.ansatz.nbetas
        
        # Split C and Emo to alpha and beta section
        Ca,Cb = np.hsplit(C,2)
        Emo_a,Emo_b = np.hsplit(Emo,2)

        # Kernel of MP2 algorithm
        results = mp2(Eri,Ca,Cb,Emo_a,Emo_b,na,nb)
        Enuc = float(ham['Nri'])
        Emp0 = results['Emp0']
        Emp1 = results['Emp1']
        Emp2 = results['Emp2']
        Etot = Enuc+Emp0+Emp1+Emp2

        log.hline()
        log('Results'.format())
        log.hline()
        log('{0:<20s} {1:<20.6f}'.format('Enuc', Enuc))
        log('{0:<20s} {1:<20.6f}'.format('Emp0', Emp0))
        log('{0:<20s} {1:<20.6f}'.format('Emp1', Emp1))
        log('{0:<20s} {1:<20.6f}'.format('Emp2', Emp2))
        log('{0:<20s} {1:<20.6f}'.format('Etot', Etot))
        log.blank()

        return Etot

    def assign_hamiltonian(self,ham):
        """Assign the chemical Hamiltonian to the solver.
        
        Parameters
        ----------
        ham: Hamiltonian
            Chemical Hamiltonian.
        """
        from moha.system.operator.hamiltonian import Hamiltonian
        if not isinstance(ham, Hamiltonian):
            raise TypeError("Parameter ham must be instance of Hamiltonian.")
        self.ham = ham

    def assign_wavefunction(self,wfn):
        """Assign the Hartree Fock wavefunction to the solver.
        
        Parameters
        ----------
        wfn: HFWaveFunction
            Hartree Fock wavefunction.
        """
        from moha.system.wavefunction.hf_wavefunction import HFWaveFunction
        if not isinstance(wfn, HFWaveFunction):
            raise TypeError("Parameter wfn must be instance of HFWaveFuntion.")
        self.wfn = wfn
