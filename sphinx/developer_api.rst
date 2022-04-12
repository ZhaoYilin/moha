API Doucmentation
#################
This detailed reference lists all the classes and functions contained in the package. If you are just looking to get started, read the :doc:`/tutorial/index` first.

The :class:`.Lattice` describes the unit cell of a crystal, while the :class:`.Model` is used to build up a larger system by translating the unit cell to fill a certain shape or symmetry. The model builds the Hamiltonian matrix by applying fields and other modifier parameters.



Molecule System
===============
.. autosummary::
    :toctree: developer_api

    moha.system.molecule
    moha.system.basis
    moha.system.operator

Hartree Fock
============
.. autosummary::
    :toctree: developer_api 

    moha.hf.scf

Post Hartree Fock
=================
.. autosummary::
    :toctree: developer_api 

    moha.posthf.ci.cisd
    moha.posthf.cc.ccsd
    moha.posthf.pt.mp

Property
========
.. autosummary::
    :toctree: developer_api 

    moha.property.population_analysis
    moha.property.multipole_moment
    moha.property.multipole_polarizability
    moha.property.excitation_energy
