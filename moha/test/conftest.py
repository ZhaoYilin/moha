import pytest
from pathlib import Path

from moha import *

cur_path = os.path.abspath(os.path.dirname(__file__))
root_path = cur_path[:cur_path.find('moha/')+len('moha/')]

@pytest.fixture
def system():
    """Build system of water molecule and sto-3g basis set.
    """
    mol_file = str(root_path)+'/moha/test/h2o.xyz'
    mol,orbs = IOSystem.from_file(mol_file,'sto-3g.nwchem')
    
    ham = ChemicalHamiltonian.build(mol,orbs)
    wfn = HFWaveFunction(10,7,{'alpha':5,'beta':5})

    return mol, orbs, ham, wfn

@pytest.fixture
def hartree_fock(system):
    """Build system of water molecule and sto-3g basis set.
    """
    mol, orbs, ham, wfn = system
    hf_results = PlainSCFSolver(ham,wfn).kernel()
    
    return hf_results
