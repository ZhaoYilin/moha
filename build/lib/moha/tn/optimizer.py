import numpy as np
import moha.tn.auxiliary as auxiliary 

class Optimizer(object):
    """
    """
    def __init__(self):
        pass

class Time_evolution(Optimizer):
    """
    Parameters:
        t_interval: tim interval for evolution
        t: how long the system evolves for
    """
    def __init__(self,Psi,H,t_interval,t,cutoff=1E-5):
        self.Psi = Psi
        self.H = H
        self.t_interval = t_interval
        self.t = t

    def tmps(self):
        pass

    def tebd(self):
        """
        Time Evolving Block Decimation(TEBD)
        """
        D = self.Psi.D
        d = self.Psi.d
        N = self.Psi.N
        #readin wavefunction in vidal form
        vertices = self.Psi.vertices
        bonds = self.Psi.bonds
        #Build the two sites time evolution operator
        TEO = auxiliary.TEO_two_sites(self.H.vertices,self.t_interval)
        # perform the iaginary time evolution 
        for step in range(0,self.t):
            for i in (range(N)[::2]+range(N)[1::2]):
                theta = auxiliary.TE_two_sites(bonds,vertices,TEO,i,N,d,D)
            print("E_iTEBD=", -np.log(np.sum(theta**2))/self.t_interval/2.)
         

class Variational(Optimizer):
    """
    """
    def __init__(self, Psi, H, m, sweeps,cutoff=1E-5):
        self.Psi = Psi
        self.H = H
        self.m = m
        self.sweeps = sweeps
        self.cutoff = cutoff

    def dmrg(self,type='two_sites'):
        """
        two-site DMRG
        """
        MPS = self.Psi.vertices
        MPO = self.H.vertices
        sweeps = self.sweeps
        m = self.m
        N = len(MPS)
        if type=='two_sites':
            EL = auxiliary.EL(MPS, MPO, MPS, 0)
            ER = auxiliary.ER(MPS, MPO, MPS, 1)
            for sweep in range(0,int(sweeps/2)):
                for i in range(0, N-2):
                    Energy,MPS[i],MPS[i+1],trunc,states = auxiliary.optimize_two_sites(MPS[i],MPS[i+1],
                                                                            MPO[i],MPO[i+1],
                                                                            EL[-1], ER[-1], m, 'left')
                    print("Sweep {:} Sites {:},{:}    Energy {:16.12f}    States {:4} Truncation {:16.12f}"
                            .format(sweep*2,i,i+1, Energy, states, trunc))
                    EL.append(auxiliary.EL_update(MPO[i], MPS[i], EL[-1], MPS[i]))
                    ER.pop();
                for i in range(N-2, 0, -1):
                    Energy,MPS[i],MPS[i+1],trunc,states = auxiliary.optimize_two_sites(MPS[i],MPS[i+1],
                                                                            MPO[i],MPO[i+1],
                                                                            EL[-1], ER[-1], m, 'right')
                    print("Sweep {} Sites {},{}    Energy {:16.12f}    States {:4} Truncation {:16.12f}"
                            .format(sweep*2+1,i,i+1, Energy, states, trunc))
                    ER.append(auxiliary.ER_update(MPO[i+1], MPS[i+1], ER[-1], MPS[i+1]))
                    EL.pop();
        elif type=='one_site':
            EL = auxiliary.EL(MPS, MPO, MPS, 0)
            ER = auxiliary.ER(MPS, MPO, MPS, 0)
            for sweep in range(0,int(sweeps/2)):
                for i in range(0, N-1):
                    Energy,MPS[i],MPS[i+1] = auxiliary.optimize_one_site(MPS[i],MPS[i+1],
                                                                                    MPO[i],EL[-1],ER[-1],
                                                                                    'left')
                    print("Sweep {:} Sites {:}    Energy {:16.12f}"
                            .format(sweep*2,i, Energy))
                    EL.append(auxiliary.EL_update(MPO[i], MPS[i], EL[-1], MPS[i]))
                    ER.pop();
                for i in range(N-1, 0, -1):
                    Energy,MPS[i],MPS[i-1] = auxiliary.optimize_one_site(MPS[i],MPS[i-1],
                                                                                    MPO[i],EL[-1], ER[-1],
                                                                                    'right')
                    print("Sweep {:} Sites {:}    Energy {:16.12f}"
                            .format(sweep*2,i, Energy))
                    ER.append(auxiliary.ER_update(MPO[i], MPS[i], ER[-1], MPS[i]))
                    EL.pop();
        self.Psi.vertices = MPS

class Projected(Optimizer):
    """
    """
    def __init__(self, Psi, H, configurations, iters, cutoff=1E-5):
        self.Psi = Psi
        self.H = H
        self.configurations = configurations
        self.iters = iters
        self.cutoff = cutoff

    def pmps(self):
        """
        pmps: projected matrix product states
        """
        MPO = self.H.vertices
        MPS = self.Psi.vertices
        for i in range(self.iters):
            MPS = auxiliary.pmps_step(self.configurations,MPO,MPS)
            print(auxiliary.Energy(MPS,MPO))

class Peturbation(Optimizer):
    """
    """
    def __init__(self):
        pass

    def pt_dmrg(self):
        pass
