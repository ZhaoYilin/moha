import moha
import numpy as np

mol,orbs = moha.IOSystem.from_file('../data/water.xyz','sto-3g.nwchem')
ham = moha.ChemicalHamiltonian.build(mol,orbs)
wfn = moha.HFWaveFunction(10,7,{'alpha':5,'beta':5})
hf_results = moha.PlainSCFSolver(ham,wfn).kernel()
