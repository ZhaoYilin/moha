from moha.system.basis.basis import SlaterDeterminant
from moha.system.wavefunction.hf_wavefunction import HFWaveFunction
from moha.io.log import log,timer

import numpy as np

__all__ = ['HFSolver']

class HFSolver(object):
    """Plain self-consistent mean field method solver.

    Attributes
    ----------
    ham: Hamiltonian
        Chemical Hamiltonian.

    sym: Symmetry
        Symmetry of the target state.

    I_thld: int
        Threshold of iteration number.

    E_thld: float
        Threshold of energy convergence.

    D_thld: float
            Threshold of density matrix convergence.

    Methods
    -------
    __init__(self,ham,sym,I_thld=50,E_thld=1.0E-8,D_thld=1.0E-8)
        Initialize the solver.

    rhf(self,acceleration=False)
        Restricted Hartree Fock solver.

    uhf(self,acceleration=False)
        Unrestricted Hartree Fock solver.
    
    assign_hamiltonian(self,ham)
        Assign the chemical Hamiltonian to the solver.

    assign_symmetry(self,sym)
        Assign symmetry to the solver.
    
    assign_iteration_threshold(self,I_thld)
        Assign the threshold of iteration number to the solver.

    assign_energy_threshold(self,E_thld)
        Assign the threshold of energy convergence to the solver.
    
    assign_density_threshold(self,D_thld)
        Assign the threshold of density convergence to the solver.
    """
    def __init__(self,ham,sym,I_thld=100,E_thld=1.0E-8,D_thld=1.0E-8):
        """Initialize the solver.

        Parameters
        ----------
        ham: Hamiltonian
            Chemical Hamiltonian.

        sym: Symmetry
            Symmetry of the target state.

        I_thld: int
            Threshold of iteration number.

        E_thld: float
            Threshold of energy convergence.

        D_thld: float
            Threshold of density matrix convergence.
        """
        self.assign_hamiltonian(ham)
        self.assign_symmetry(sym)
        self.assign_iteration_threshold(I_thld)
        self.assign_energy_threshold(E_thld)
        self.assign_density_threshold(D_thld)    
    
    @timer.with_section('HF')
    def hf(self,acceleration=False):
        """Hartree Fock solver.

        Parameters
        ----------
        acceleration : bool
            Acceleration argument.

        Returns
        -------
        E : float
            Hartree Fock energy.

        wfn : HFWaveFunction
            Hartree Fock wavefcuntion.
        """
        if self.symmetry.ms2==0:
            E,wfn = self.rhf(acceleration)
        else:
            E,wfn = self.uhf(acceleration)

        return E,wfn
        


    @timer.with_section('HF')
    def rhf(self,acceleration=False):
        """Restricted Hartree Fock solver.

        Parameters
        ----------
        acceleration : bool
            Acceleration argument.

        Returns
        -------
        E : float
            Hartree Fock energy.

        wfn : HFWaveFunction
            Hartree Fock wavefcuntion.
        """
        log.hline(char='=')
        log.blank()
        log('Hartree Fock Section'.format())

        #Set attributes
        Nri = float(self.ham['Nri'])
        S = self.ham['S']
        Hcore = self.ham['Hcore']
        Eri = self.ham['Eri']

        norb = S.shape[0]
        ne = self.symmetry.n
        ndocc = int(ne/2)
        ms2 = self.symmetry.ms2
        if not (ne%2==0 and ms2==0):
            raise ValueError("System must be a restricted close shell.")        

        I_thld = self.I_thld
        E_thld = self.E_thld
        D_thld = self.D_thld
        
        if acceleration:
            from moha.hf.rhf_diis import rhf_diis
            results = rhf_diis(S,Hcore,Eri,ndocc,I_thld,E_thld,D_thld)
        else:
            from moha.hf.rhf import rhf
            results = rhf(S,Hcore,Eri,ndocc,I_thld,E_thld,D_thld)

        log.hline()
        log('Results'.format())
        log.hline()
        log('{0:<20s} {1:<20.6f}'.format('Enuc', Nri))
        log('{0:<20s} {1:<20.6f}'.format('Eele', results['energy']))
        log('{0:<20s} {1:<20.6f}'.format('Etot', results['energy']+Nri))
        log.blank()

        E = results['energy']

        wfn = HFWaveFunction()
        ansatz = SlaterDeterminant.ground(ne,2*norb)
        wfn.assign_ansatz(ansatz)
        wfn.assign_coefficient(results['C'])
        wfn.assign_density_matrix(results['D'])
        wfn.assign_orbital_energies(results['Emo'])
        
        return E,wfn

    @timer.with_section('HF')
    def uhf(self,acceleration=False):
        """Unrestricted Hartree Fock solver.

        Parameters
        ----------
        acceleration : bool
            Acceleration argument.

        Returns
        -------
        E : float
            Hartree Fock energy.

        wfn : HFWaveFunction
            Hartree Fock wavefcuntion.
        """
        log.hline(char='=')
        log.blank()
        log('Hartree Fock Section'.format())

        #Set attributes
        Nri = float(self.ham['Nri'])
        S = self.ham['S']
        Hcore = self.ham['Hcore']
        Eri = self.ham['Eri']

        norb = S.shape[0]
        ne = self.symmetry.n
        ms2 = self.symmetry.ms2
        na = int((ne+ms2)/2)
        nb = int((ne-ms2)/2)

        I_thld = self.I_thld
        E_thld = self.E_thld
        D_thld = self.D_thld
        
        if acceleration:
            from moha.hf.uhf_diis import uhf_diis
            results = uhf_diis(S,Hcore,Eri,na,nb,I_thld,E_thld,D_thld)
        else:
            from moha.hf.uhf import uhf
            results = uhf(S,Hcore,Eri,na,nb,I_thld,E_thld,D_thld)

        log.hline()
        log('Results'.format())
        log.hline()
        log('{0:<20s} {1:<20.6f}'.format('Enuc', Nri))
        log('{0:<20s} {1:<20.6f}'.format('Eele', results['energy']))
        log('{0:<20s} {1:<20.6f}'.format('Etot', results['energy']+Nri))
        log.blank()

        E = results['energy']

        wfn = HFWaveFunction()
        alpha = [1]*na+[0]*(norb-na)
        beta = [1]*nb+[0]*(norb-nb)
        configuration = []
        for i,j in zip(alpha, beta):
            configuration.append(i)
            configuration.append(j)
        ansatz = SlaterDeterminant(configuration)
        wfn.assign_ansatz(ansatz)
        wfn.assign_coefficient(results['C'])
        wfn.assign_density_matrix(results['D'])
        wfn.assign_orbital_energies(results['Emo'])
        
        return E,wfn

    def assign_hamiltonian(self,ham):
        """Assign the chemical Hamiltonian to the solver.
        
        Parameters
        ----------
        ham: Hamiltonian
            Chemical Hamiltonian.

        Raises
        ------
        TypeError
            If ham is not Hamiltonian instance.
        """
        from moha.system.operator.hamiltonian import Hamiltonian
        if not isinstance(ham, Hamiltonian):
            raise TypeError("Parameter ham must be Hamiltonian instance.")
        self.ham = ham

    def assign_symmetry(self,sym):
        """Assign symmetry to the solver.

        Parameters
        ----------
        sym: Symmetry
            Symmetry of the target state.

        Raises
        ------
        TypeError
            If symmetry is not Symmetry instance.
        """
        from moha.symmetry.symmetry import Symmetry
        if not isinstance(sym,Symmetry):
            raise TypeError("Parameters sym must be Symmetry instance.")
        self.symmetry = sym

    def assign_iteration_threshold(self,I_thld):
        """Assign the threshold of iteration number to the solver.

        Parameters
        ----------
        I_thld: int
            Threshold of iteration number.
        
        Raises
        ------
        TypeError
            If threshold of iteration number is not integer.

        ValueError
            If threshold of iteration number is not positive.
        """
        if not isinstance(I_thld, int):
            raise TypeError("Threshold of iteration number must be integer.")
        if I_thld <= 0:
            raise ValueError("Threshold of iteration number must be positive.")
        self.I_thld = I_thld       

    def assign_energy_threshold(self,E_thld):
        """Assign the threshold of energy convergence to the solver.

        Parameters
        ----------
        E_thld: float
            Threshold of energy convergence.
        
        Raises
        ------
        TypeError
            If threshold of energy convergence is not integer.

        ValueError
            If threshold of energy convergence is not positive.
        """
        if not isinstance(E_thld, float):
            raise TypeError("Threshold of energy convergence must be float.")
        if E_thld <= 0:
            raise ValueError("Threshold of energy convergence must be positive.")        
        self.E_thld = E_thld

    def assign_density_threshold(self,D_thld):
        """Assign the threshold of density convergence to the solver.

        Parameters
        ----------
        D_thld: float
            Threshold of density convergence.
        
        Raises
        ------
        TypeError
            If threshold of density convergence is not integer.

        ValueError
            If threshold of density convergence is not positive.
        """
        if not isinstance(D_thld, float):
            raise TypeError("Threshold of density convergence must be float.")
        if D_thld <= 0:
            raise ValueError("Threshold of density convergence must be positive.")        
        self.D_thld = D_thld
