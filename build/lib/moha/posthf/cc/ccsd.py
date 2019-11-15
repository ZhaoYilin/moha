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

def initialize(fs,spinints,Nelec,dim):
    """
    Init empty T1 (ts) and T2 (td) arrays
    and make denominator arrays Dai, Dabij
    """
    # Initial guess for T1 and T2
    ts = np.zeros((dim,dim))
    td = np.zeros((dim,dim,dim,dim))
    for a in range(Nelec,dim):
        for b in range(Nelec,dim):
            for i in range(0,Nelec):
                for j in range(0,Nelec):
                    td[a,b,i,j] += spinints[i,j,a,b]/(fs[i,i] + fs[j,j] - fs[a,a] - fs[b,b])

    # Equation (12) of Stanton
    Dai = np.zeros((dim,dim))
    for a in range(Nelec,dim):
        for i in range(0,Nelec):
            Dai[a,i] = fs[i,i] - fs[a,a]

    # Stanton eq (13)
    Dabij = np.zeros((dim,dim,dim,dim))
    for a in range(Nelec,dim):
        for b in range(Nelec,dim):
            for i in range(0,Nelec):
                for j in range(0,Nelec):
                    Dabij[a,b,i,j] = fs[i,i] + fs[j,j] - fs[a,a] - fs[b,b]
    return ts,td,Dai,Dabij

# Stanton eq (9)
def taus(ts,td,a,b,i,j):
  taus = td[a,b,i,j] + 0.5*(ts[a,i]*ts[b,j] - ts[b,i]*ts[a,j])
  return taus

# Stanton eq (10)
def tau(ts,td,a,b,i,j):
  tau = td[a,b,i,j] + ts[a,i]*ts[b,j] - ts[b,i]*ts[a,j]
  return tau

# We need to update our intermediates at the beginning, and 
# at the end of each iteration. Each iteration provides a new
# guess at the amplitudes T1 (ts) and T2 (td), that *hopefully*
# converges to a stable, ground-state, solution.

def updateintermediates(x,Nelec,dim,fs,spinints,ts,td):
  if x == True:
    # Stanton eq (3)
    Fae = np.zeros((dim,dim))
    for a in range(Nelec,dim):
      for e in range(Nelec,dim):
        Fae[a,e] = (1 - (a == e))*fs[a,e]
        for m in range(0,Nelec):
          Fae[a,e] += -0.5*fs[m,e]*ts[a,m]
          for f in range(Nelec,dim):
            Fae[a,e] += ts[f,m]*spinints[m,a,f,e] 
            for n in range(0,Nelec):
              Fae[a,e] += -0.5*taus(ts,td,a,f,m,n)*spinints[m,n,e,f]

    # Stanton eq (4)
    Fmi = np.zeros((dim,dim))
    for m in range(0,Nelec):
      for i in range(0,Nelec):
        Fmi[m,i] = (1 - (m == i))*fs[m,i]
        for e in range(Nelec,dim):
          Fmi[m,i] += 0.5*ts[e,i]*fs[m,e]
          for n in range(0,Nelec):
            Fmi[m,i] += ts[e,n]*spinints[m,n,i,e] 
            for f in range(Nelec,dim):
              Fmi[m,i] += 0.5*taus(ts,td,e,f,i,n)*spinints[m,n,e,f]

    # Stanton eq (5)
    Fme = np.zeros((dim,dim))
    for m in range(0,Nelec):
      for e in range(Nelec,dim):
        Fme[m,e] = fs[m,e]
        for n in range(0,Nelec):
          for f in range(Nelec,dim):
            Fme[m,e] += ts[f,n]*spinints[m,n,e,f]

    # Stanton eq (6)
    Wmnij = np.zeros((dim,dim,dim,dim))
    for m in range(0,Nelec):
      for n in range(0,Nelec):
        for i in range(0,Nelec):
          for j in range(0,Nelec):
            Wmnij[m,n,i,j] = spinints[m,n,i,j]
            for e in range(Nelec,dim):
              Wmnij[m,n,i,j] += ts[e,j]*spinints[m,n,i,e] - ts[e,i]*spinints[m,n,j,e]
              for f in range(Nelec,dim):
                Wmnij[m,n,i,j] += 0.25*tau(ts,td,e,f,i,j)*spinints[m,n,e,f]

    # Stanton eq (7)
    Wabef = np.zeros((dim,dim,dim,dim))
    for a in range(Nelec,dim):
      for b in range(Nelec,dim):
        for e in range(Nelec,dim):
          for f in range(Nelec,dim):
            Wabef[a,b,e,f] = spinints[a,b,e,f]
            for m in range(0,Nelec):
              Wabef[a,b,e,f] += -ts[b,m]*spinints[a,m,e,f] + ts[a,m]*spinints[b,m,e,f]
              for n in range(0,Nelec):
                Wabef[a,b,e,f] += 0.25*tau(ts,td,a,b,m,n)*spinints[m,n,e,f]

    # Stanton eq (8)
    Wmbej = np.zeros((dim,dim,dim,dim))
    for m in range(0,Nelec):
      for b in range(Nelec,dim):
        for e in range(Nelec,dim):
          for j in range(0,Nelec):
            Wmbej[m,b,e,j] = spinints[m,b,e,j]
            for f in range(Nelec,dim):
              Wmbej[m,b,e,j] += ts[f,j]*spinints[m,b,e,f]
            for n in range(0,Nelec):
              Wmbej[m,b,e,j] += -ts[b,n]*spinints[m,n,e,j]
              for f in range(Nelec,dim):
                Wmbej[m,b,e,j] += -(0.5*td[f,b,j,n] + ts[f,j]*ts[b,n])*spinints[m,n,e,f]

    return Fae, Fmi, Fme, Wmnij, Wabef, Wmbej

# makeT1 and makeT2, as they imply, construct the actual amplitudes necessary for computing
# the CCSD energy (or computing an EOM-CCSD Hamiltonian, etc)

# Stanton eq (1)
def makeT1(x,Nelec,dim,fs,spinints,ts,td,Dai,Fae,Fmi,Fme):
  if x == True:
    tsnew = np.zeros((dim,dim))
    for a in range(Nelec,dim):
      for i in range(0,Nelec):
        tsnew[a,i] = fs[i,a]
        for e in range(Nelec,dim):
          tsnew[a,i] += ts[e,i]*Fae[a,e]
        for m in range(0,Nelec):
          tsnew[a,i] += -ts[a,m]*Fmi[m,i]
          for e in range(Nelec,dim):
            tsnew[a,i] += td[a,e,i,m]*Fme[m,e]
            for f in range(Nelec,dim):
              tsnew[a,i] += -0.5*td[e,f,i,m]*spinints[m,a,e,f]
            for n in range(0,Nelec):
              tsnew[a,i] += -0.5*td[a,e,m,n]*spinints[n,m,e,i]
        for n in range(0,Nelec):
          for f in range(Nelec,dim): 
            tsnew[a,i] += -ts[f,n]*spinints[n,a,i,f]
        tsnew[a,i] = tsnew[a,i]/Dai[a,i]
  return tsnew

# Stanton eq (2)
def makeT2(x,Nelec,dim,fs,spinints,ts,td,Dabij,Fae,Fmi,Fme,Wmnij,Wabef,Wmbej):
  if x == True:
    tdnew = np.zeros((dim,dim,dim,dim))
    for a in range(Nelec,dim):
      for b in range(Nelec,dim):
        for i in range(0,Nelec):
          for j in range(0,Nelec):
            tdnew[a,b,i,j] += spinints[i,j,a,b]
            for e in range(Nelec,dim):
              tdnew[a,b,i,j] += td[a,e,i,j]*Fae[b,e] - td[b,e,i,j]*Fae[a,e]
              for m in range(0,Nelec):
                tdnew[a,b,i,j] += -0.5*td[a,e,i,j]*ts[b,m]*Fme[m,e] + 0.5*td[a,e,i,j]*ts[a,m]*Fme[m,e]
                continue
            for m in range(0,Nelec):
              tdnew[a,b,i,j] += -td[a,b,i,m]*Fmi[m,j] + td[a,b,j,m]*Fmi[m,i]
              for e in range(Nelec,dim):
                tdnew[a,b,i,j] += -0.5*td[a,b,i,m]*ts[e,j]*Fme[m,e] + 0.5*td[a,b,i,m]*ts[e,i]*Fme[m,e]
                continue
            for e in range(Nelec,dim):
              tdnew[a,b,i,j] += ts[e,i]*spinints[a,b,e,j] - ts[e,j]*spinints[a,b,e,i]
              for f in range(Nelec,dim):
                tdnew[a,b,i,j] += 0.5*tau(ts,td,e,f,i,j)*Wabef[a,b,e,f]
                continue
            for m in range(0,Nelec):
              tdnew[a,b,i,j] += -ts[a,m]*spinints[m,b,i,j] + ts[b,m]*spinints[m,a,i,j]  
              for e in range(Nelec,dim):
                tdnew[a,b,i,j] += td[a,e,i,m]*Wmbej[m,b,e,j] - ts[e,i]*ts[a,m]*spinints[m,b,e,j]
                tdnew[a,b,i,j] += -td[a,e,j,m]*Wmbej[m,b,e,i] + ts[e,j]*ts[a,m]*spinints[m,b,e,i]
                tdnew[a,b,i,j] += -td[b,e,i,m]*Wmbej[m,a,e,j] - ts[e,i]*ts[b,m]*spinints[m,a,e,j]
                tdnew[a,b,i,j] += td[b,e,j,m]*Wmbej[m,a,e,i] - ts[e,j]*ts[b,m]*spinints[m,a,e,i]
                continue
              for n in range(0,Nelec):
                tdnew[a,b,i,j] += 0.5*tau(ts,td,a,b,m,n)*Wmnij[m,n,i,j]
                continue
            tdnew[a,b,i,j] = tdnew[a,b,i,j]/Dabij[a,b,i,j] 
    return tdnew

# Expression from Crawford, Schaefer (2000) 
# DOI: 10.1002/9780470125915.ch2
# Equation (134) and (173)
# computes CCSD energy given T1 and T2
def ccsdenergy(Nelec,dim,fs,spinints,ts,td):
  ECCSD = 0.0
  for i in range(0,Nelec):
    for a in range(Nelec,dim):
      ECCSD += fs[i,a]*ts[a,i]
      for j in range(0,Nelec):
        for b in range(Nelec,dim):
          ECCSD += 0.25*spinints[i,j,a,b]*td[a,b,i,j] + 0.5*spinints[i,j,a,b]*(ts[a,i])*(ts[b,j]) 
  return ECCSD
#######################################################
#
#   CCSD CALCULATION
#
#######################################################

class CCSolver(object):
    def __init__(self,maxiter,cutoff):
        self.maxiter = maxiter
        self.cutoff = cutoff

    @classmethod
    def ccsd(cls,hfwavefunction,hamiltonian,maxiter=100,cutoff=1e-6):
        occ = hfwavefunction.occ
        Nelec = occ['alpha'] + occ['beta']
        C = hfwavefunction.coefficient
        dim = hamiltonian.dim*2
        #Transfer Fock integral from spatial to spin basis
        fs = spinfock(hfwavefunction.eorbitals)
        #Transfer electron repulsion integral from atomic basis
        #to molecular basis
        hamiltonian.operators['electron_repulsion'].basis_transformation(C)
        #build double bar integral <ij||kl>
        spinints = hamiltonian.operators['electron_repulsion'].double_bar 
        ts,td,Dai,Dabij = initialize(fs,spinints,Nelec,dim)
        ECCSD = 0
        DECC = 1.0
        Iter = 0
        log.hline()
        log('{0:2s} {1:3s} {2:4s}'.format('Iter', 'ECC', 'Delta'))
        log.hline()
        while DECC > cutoff or Iter> maxiter: # arbitrary convergence criteria
            Iter += 1
            OLDCC = ECCSD
            Fae,Fmi,Fme,Wmnij,Wabef,Wmbej = updateintermediates(True,Nelec,dim,fs,spinints,ts,td)
            ts = makeT1(True,Nelec,dim,fs,spinints,ts,td,Dai,Fae,Fmi,Fme)
            td = makeT2(True,Nelec,dim,fs,spinints,ts,td,Dabij,Fae,Fmi,Fme,Wmnij,Wabef,Wmbej)
            ECCSD = ccsdenergy(Nelec,dim,fs,spinints,ts,td)
            DECC = abs(ECCSD - OLDCC)
            log('{0:2d} {1:3f} {2:4f}'.format(Iter, ECCSD, DECC))
        print('CC iterations converged')
        print("E(corr,CCSD) = ", ECCSD)
