from moha import *
mol,orbs = IOSystem.from_file('h2o.xyz','sto-3g.nwchem')
ham = RestrictedChemicalHamiltonian.build(mol,orbs)
wfn = HFWaveFunction(10,7,{'alpha':5,'beta':5})
hf_results = DIISSCFSolver(ham,wfn).kernel()
print(hf_results['orbital_energies'])
