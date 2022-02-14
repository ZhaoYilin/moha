from moha import *

ovlp = np.load('./data/overlap.npy')
h1e = np.load('./data/kinetic.npy')+np.load('./data/nuclear_attraction.npy')
g2e = np.load('./data/electron_repulsion.npy')

ham = ChemicalHamiltonian.from_numpy(8.002367,ovlp,h1e,g2e)

wfn = HFWaveFunction(10,7,{'alpha':5,'beta':5})

scfsolver = PlainSCFSolver(ham,wfn)
hf_results = scfsolver.kernel()

