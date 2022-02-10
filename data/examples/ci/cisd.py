from moha import *

mol,orbs = IOSystem.from_file('../../molecules/water.xyz','sto-3g.nwchem')
ham = ChemicalHamiltonian.build(mol,orbs)
wfn = HFWaveFunction(10,7,{'alpha':5,'beta':5})

scfsolver = PlainSCFSolver(ham,wfn)
hf_results = scfsolver.kernel()

cisolver = CISDSolver(ham,wfn,hf_results)
ci_results = cisolver.kernel()
