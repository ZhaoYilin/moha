from moha.modelsystem import sites
from moha.tn import tns
from moha.tn import tnso
from moha.tn import optimizer
import numpy as np
Site = sites.SpinOneHalfSite()
### initial state |+-+-+-+-+->
Psi = tns.MPS(Site,24,4,'MPS')
### Hamiltonian MPO5
H = tnso.MPO(Site,24,'Heisenberg',{'J':1.0,'Jz':1.0,'h':0.0})
Optimizer = optimizer.Variational(Psi,H,50,100)
Optimizer.dmrg()
