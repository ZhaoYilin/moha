from moha import *
mol,orbs = IOSystem.from_file('h2o.xyz','sto-3g.nwchem')
ham = Hamiltonian.build(mol,orbs)
hfwavefunction = SCFSolver.hf(ham,occ = {'alpha':5,'beta':5})
PopulationAnalysis.mulliken(mol,orbs,hfwavefunction,ham)
