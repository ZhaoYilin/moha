from moha.property.auxiliary import *
from moha.io.log import log, timer

import numpy as np

__all__ = ['ExcitationEnergyCIS','ExcitationEnergyTDHF']

class ExcitationEnergy(object):
    """Excitation energy solver class.

    Attributes
    ----------
    ham
        Chemical Hamiltonian.

    wfn
        Hartree Fock wavefunction.
    
    level : {None, int, list}
        Excitation energy level.

    Methods
    -------
    __init__(self,ham,wfn)
        Initialize the solver.

    assign_hamiltonian(self,ham)
        Assign the chemical Hamiltonian to the solver.

    assign_wavefunction(self,wfn)
        Assign the Hartree Fock wavefunction to the solver.
    
    assign_level(self,level)
        Assign the excitation energy level to the solver.
    """
    def __init__(self,ham,wfn):
        """Initialize the solver.

        Attributes
        ----------
        ham
            Chemical Hamiltonian.

        wfn
            Hartree Fock wavefunction.
        """
        self.assign_hamiltonian(ham)
        self.assign_wavefunction(wfn)

    def assign_hamiltonian(self,ham):
        """Assign the chemical Hamiltonian to the solver.

        Attributes
        ----------
        ham
            Chemical Hamiltonian.
        """
        self.ham = ham

    def assign_wavefunction(self,wfn):
        """Assign the Hartree Fock wavefunction to the solver.

        Attributes
        ----------
        wfn
            Hartree Fock wavefunction.
        """
        self.wfn = wfn    

    def assign_level(self,level):
        """Assign the excitation energy level to the solver.

        Attributes
        ----------
        level
            Excitation energy level.

        Raises
        ------
        TypeError
            If excitation energy level is not an integer.
        ValueError
            If excitation energy level is not a positive number.

        """
        if level is None:
            self.level = level
        if not isinstance(level, (int,list)):
            raise TypeError("Excitation energy level must be int or list")
        self.level = level-1    

class ExcitationEnergyCIS(ExcitationEnergy):
    """Excitation energy solver using configuration interaction single method.

    Attributes
    ----------
    ham
        Chemical Hamiltonian.

    wfn
        Hartree Fock wavefunction.
    
    level : {None, int, list}
        Excitation energy level.

    Methods
    -------
    __init__(self,ham,wfn)
        Initialize the solver.

    kernel(self,level=None)
        Kernel of the solver.

    assign_hamiltonian(self,ham)
        Assign the chemical Hamiltonian to the solver.

    assign_wavefunction(self,wfn)
        Assign the Hartree Fock wavefunction to the solver.
    
    assign_level(self,level)
        Assign the excitation energy level to the solver.
    """
    def __init__(self,ham,wfn):
        """Initialize the solver.

        Attributes
        ----------
        ham
            Chemical Hamiltonian.

        wfn
            Hartree Fock wavefunction.
        """
        super().__init__(ham,wfn)

    @timer.with_section('Excitation')
    def kernel(self,level=None):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            Excitation energy calculation results.
        """    
        log.hline()
        log('Excitation Energy Calculation Section'.format())
        log.hline()

        occ = self.wfn.occ
        Nelec = occ['alpha'] + occ['beta']
        dim = self.ham.nspatial
        C = self.wfn.coefficients
        Nov = Nelec*(2*dim - Nelec)
        #Transfer Fock integral from spatial to spin basis
        fs = spinfock(self.wfn.orbital_energies)
        #Transfer electron repulsion integral from atomic basis
        #to molecular basis
        self.ham.operators['electron_repulsion'].basis_transformation(C)
        #build double bar integral <ij||kl>
        spinints = self.ham.operators['electron_repulsion'].double_bar
        
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
        ECIS = np.sort(ECIS)[level].flatten()

        log.hline()
        log('Excitaiton Energy by CIS'.format())
        for i,energy in enumerate(ECIS):
            log('{0:2d}\t{1:6f}'.format(i,energy))
        log.hline()

        results = {
        "success": True,
        "excitation_energies": ECIS
        }
        return results

class ExcitationEnergyTDHF(ExcitationEnergy):
    """Excitation energy solver using time dependent Hartree Fock method.

    Attributes
    ----------
    ham
        Chemical Hamiltonian.

    wfn
        Hartree Fock wavefunction.
    
    level : {None, int, list}
        Excitation energy level.

    Methods
    -------
    __init__(self,ham,wfn)
        Initialize the solver.

    kernel(self,level=None)
        Kernel of the solver.

    assign_hamiltonian(self,ham)
        Assign the chemical Hamiltonian to the solver.

    assign_wavefunction(self,wfn)
        Assign the Hartree Fock wavefunction to the solver.
    
    assign_level(self,level)
        Assign the excitation energy level to the solver.
    """
    def __init__(self,ham,wfn):
        """Initialize the solver.

        Attributes
        ----------
        ham
            Chemical Hamiltonian.

        wfn
            Hartree Fock wavefunction.
        """
        super().__init__(ham,wfn)

    @timer.with_section('Excitation')
    def kernel(self,level=None):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            Excitation energy calculation results.
        """    
        

        log.hline()
        log('Excitation Energy Calculation Section'.format())
        log.hline()

        occ = self.wfn.occ
        Nelec = occ['alpha'] + occ['beta']
        dim = self.ham.nspatial
        C = self.wfn.coefficients
        Nov = Nelec*(2*dim - Nelec)
        #Transfer Fock integral from spatial to spin basis
        fs = spinfock(self.wfn.orbital_energies)
        #Transfer electron repulsion integral from atomic basis
        #to molecular basis
        self.ham.operators['electron_repulsion'].basis_transformation(C)
        #build double bar integral <ij||kl>
        spinints = self.ham.operators['electron_repulsion'].double_bar

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
        ETD = np.sort(ETD)[level].flatten()
        
        log.hline()
        log('Excitaiton Energy by TDHF'.format())
        log('{}'.format(ETD))
        log.hline()

        results = {
        "success": True,
        "excitation_energies": ETD
        }
        return results

