from moha import *

geo = [[8,   0.000000000000,  -0.143225816552,   0.000000000000],
    [1,   1.638036840407,   1.136548822547,  -0.000000000000],
    [1,  -1.638036840407,   1.136548822547,  -0.000000000000]]

mol = Molecule.build(geo)
orb = BasisSet.build(mol,'sto-3g.nwchem')

ham = ChemicalHamiltonian.build(mol,orb)
wfn = HFWaveFunction(10,7,{'alpha':5,'beta':5})

hf_solver = PlainSCFSolver(ham,wfn)
hf_results = hf_solver.kernel()

ptsolver = MP3Solver(ham,wfn,hf_results)
pt_results = ptsolver.kernel()
