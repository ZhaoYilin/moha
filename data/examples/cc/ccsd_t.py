from moha import *

mol,orbs = IOSystem.from_file('../../molecules/water.xyz','sto-3g.nwchem')
ham = Hamiltonian.build(mol,orbs)

scfsolver = SCF_Solver()
hfw = scfsolver(ham,occ = {'alpha':5,'beta':5})

solver = CCSD_T_Solver()
solver(hfw,ham)
