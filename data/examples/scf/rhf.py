from moha import *
mol,orbs = IOSystem.from_file('water.xyz','sto-3g.nwchem')
ham = RestrictedChemicalHamiltonian.build(mol,orbs)
wfn = HFWaveFunction(10,7,{'alpha':5,'beta':5})
hf_results = PlainSCFSolver(ham,wfn).kernel()
