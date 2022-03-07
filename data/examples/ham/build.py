from moha import *
import numpy as np

mol,orbs = IOSystem.from_file('../data/water.xyz','sto-3g.nwchem')
ham = ChemicalHamiltonian.build(mol,orbs)

wfn = HFWaveFunction(10,7,{'alpha':5,'beta':5})

scfsolver = PlainSCFSolver(ham,wfn)
hf_results = scfsolver.kernel()
