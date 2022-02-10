from moha import *

mol,orbs = IOSystem.from_file('../../molecules/water.xyz','sto-3g.nwchem')
ham = ChemicalHamiltonian.build(mol,orbs)
wfn = HFWaveFunction(10,7,{'alpha':5,'beta':5})

scfsolver = PlainSCFSolver(ham,wfn)
hf_results = scfsolver.kernel()

ccsolver = CCSDSolver(ham,wfn,hf_results)
cc_results = ccsolver.kernel()
