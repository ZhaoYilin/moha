from moha import *
mol,orbs = IOSystem.from_file('h2o.xyz','sto-3g.nwchem')
ham = Hamiltonian.build(mol,orbs)
hfwavefunction = SCFSolver.hf(ham,occ = {'alpha':4,'beta':5})
ExcitationEnergy.cis(hfwavefunction,ham)
