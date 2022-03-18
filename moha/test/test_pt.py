import numpy as np
import pytest
from moha import *


def test_mp2(system, hartree_fock):
    """Test cases with mp2 slover.
    """
    mol, orbs, ham, wfn = system
    hf_results = hartree_fock
    
    ptsolver = MP2Solver(ham,wfn,hf_results)
    pt_results = ptsolver.kernel()
    mp2_energy = pt_results['mp2_energy']
    energy = pt_results['total_energy']

    assert np.allclose(mp2_energy, -0.049150)
    assert np.allclose(energy, -74.991230)

def test_mp3(system, hartree_fock):
    """Test cases with mp3 slover.
    """
    mol, orbs, ham, wfn = system
    hf_results = hartree_fock
    
    ptsolver = MP3Solver(ham,wfn,hf_results)
    pt_results = ptsolver.kernel()
    mp2_energy = pt_results['mp2_energy']
    mp3_energy = pt_results['mp3_energy']
    energy = pt_results['total_energy']

    assert np.allclose(mp2_energy, -0.049150)
    assert np.allclose(mp3_energy, -0.0002316)
    assert np.allclose(energy, -74.991461)

