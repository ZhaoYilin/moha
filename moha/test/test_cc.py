import numpy as np
import pytest
from moha import *


def test_ccsd(system, hartree_fock):
    """Test cases with ccsd slover.
    """
    mol, orbs, ham, wfn = system
    hf_results = hartree_fock
    
    ccsolver = CCSDSolver(ham,wfn,hf_results)
    cc_results = ccsolver.kernel()
    cc_energy = cc_results['total_energy']

    assert np.allclose(cc_energy, -75.012779)

def test_ccsd_t(system, hartree_fock):
    """Test cases with ccsd_t slover.
    """
    mol, orbs, ham, wfn = system
    hf_results = hartree_fock
    
    ccsolver = CCSD_TSolver(ham,wfn,hf_results)
    cc_results = ccsolver.kernel()
    cc_energy = cc_results['total_energy']

    assert np.allclose(cc_energy, -75.012879)
