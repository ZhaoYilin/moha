from moha import *
import pytest

@pytest.fixture
def system():
    """Build system of water molecule and sto-3g basis set.
    """
    geo = ([8,   0.000000000000,  -0.143225816552,   0.000000000000],
        (1,   1.638036840407,   1.136548822547,  -0.000000000000),
        [1,  -1.638036840407,   1.136548822547,  -0.000000000000])

    mol = Molecule.build(geo)
    orb = BasisSet.build(mol,'sto-3g.nwchem')
    ham = ChemicalHamiltonian.build(mol,orb)
    wfn = HFWaveFunction(10,7,{'alpha':5,'beta':5})

    return mol, orbs, ham, wfn

@pytest.fixture
def hartree_fock(system):
    """Build system of water molecule and sto-3g basis set.
    """
    mol, orbs, ham, wfn = system
    hf_results = PlainSCFSolver(ham,wfn).kernel()
    
    return hf_results
