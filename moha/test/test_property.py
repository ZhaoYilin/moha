import numpy as np
import pytest
from moha import *

def test_excitation_energy(system, hartree_fock):
    """Excitation energy test cases.
    """
    mol, orbs, ham, wfn = system
    hf_results = hartree_fock
    
    ee_cis_results = ExcitationEnergyCIS(ham,wfn).kernel()
    ee_cis_0 = ee_cis_results['excitation_energies'][0]

    assert np.allclose(ee_cis_0, 0.287256)



def test_multipolemoment(system, hartree_fock):
    """Test cases with ccsd_t slover.
    """
    mol, orbs, ham, wfn = system
    hf_results = hartree_fock

    mm_results = MultipoleMoment(mol,orbs,wfn).kernel([[1,0,0],[0,1,0],[0,0,1]])
    mms = mm_results['mms']
    
    assert np.allclose(mms, [0.000000, -4.962612, 0.000000])


def test_population(system, hartree_fock):
    """Test cases with ccsd_t slover.
    """
    mol, orbs, ham, wfn = system
    hf_results = hartree_fock

    m_results = PopulationAnalysisMulliken(mol,orbs,ham,wfn).kernel()
    l_results = PopulationAnalysisLowdin(mol,orbs,ham,wfn).kernel()
    m_population = m_results['population'] 
    l_population = l_results['population'] 
    
    assert np.allclose(m_population, [-0.253146, 0.126573, 0.126573])
    assert np.allclose(l_population, [-0.184234, 0.092117, 0.092117])
