import numpy as np
import pytest
from moha import *

def test_cis(system, hartree_fock):
    """Test cases with cis slover.
    """
    mol, orbs, ham, wfn = system
    hf_results = hartree_fock
    
    cisolver = CISSolver(ham,wfn,hf_results)
    ci_results, ciwfn = cisolver.kernel()
    ci_energy = ci_results['total_energy']

    assert np.allclose(ci_energy, -74.942080)

def test_cisd(system, hartree_fock):
    """Test cases with cisd slover.
    """
    mol, orbs, ham, wfn = system
    hf_results = hartree_fock

    cisolver = CISDSolver(ham,wfn,hf_results)
    ci_results, ciwfn = cisolver.kernel()
    ci_energy = ci_results['total_energy']

    assert np.allclose(ci_energy, -74.988112)

def test_fullci(system, hartree_fock):
    """Test cases with full ci slover.
    """
    mol, orbs, ham, wfn = system
    hf_results = hartree_fock
    
    cisolver = FullCISolver(ham,wfn,hf_results)
    ci_results, ciwfn = cisolver.kernel()
    ci_energy = ci_results['total_energy']

    assert np.allclose(ci_energy, -74.988611)
