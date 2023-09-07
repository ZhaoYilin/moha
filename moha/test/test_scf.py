from moha import *
import pytest
import numpy as np
import os

def test_plain_scf(system):
    """Test cases with plain scf slover.
    """
    mol, orbs, ham, wfn = system
    
    hf_solver = PlainSCFSolver(ham,wfn)
    hf_results = hf_solver.kernel()
    hf_energy = hf_results['total_energy']
    
    assert np.allclose(hf_energy, -74.942080)

def test_diis_scf(system):
    """Test cases with diis scf slover.
    """
    mol, orbs, ham, wfn = system
    
    hf_solver = DIISSCFSolver(ham,wfn)
    hf_results = hf_solver.kernel()
    hf_energy = hf_results['total_energy']
    
    assert np.allclose(hf_energy, -74.942080)
