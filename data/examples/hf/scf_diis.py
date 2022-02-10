from moha import *

mol,orbs = IOSystem.from_file('../../molecules/water.xyz','sto-3g.nwchem')
ham = ChemicalHamiltonian.build(mol,orbs)
wfn = HFWaveFunction(10,7,{'alpha':5,'beta':5})

hf_solver = DIISSCFSolver(ham,wfn)
hf_results = hf_solver.kernel()
