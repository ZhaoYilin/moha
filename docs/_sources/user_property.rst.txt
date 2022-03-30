Property
########

This section containts information about atomic and molecular properties that can be calculated.

Population Analysis/Atomic Charges
==================================
Population analysis is the study of charge distribution in molecules, it requires the overlap integrals and the electron density, in addition to information about the number of basis functions centered on each atom. The charge on atom A may be computed as:

Mulliken Population Analyasis
-------------------------------

.. math::
    q_A = Z_A - 2 \sum_{\mu \in A} ({\mathbf D S})_{\mu\mu},

where the summation is limited to only those basis functions centered on atom A.

.. literalinclude:: ../data/examples/property/population_analysis.py
        :lines: 1-20
        :caption: /data/examples/population_analysis.py


Multipole Moment
================
The calculation of one-electron properties requires density matrix and the
relevant property integrals. The electronic contribution to the
electric-miltipole moment may be computed using:

.. math::
    {\mu}^{klm}=-\sum_{i}\langle \Psi_0|x_i^ky_i^lz_i^m|\Psi_0 \rangle+\sum_A Z_Ax_A^ky_A^lz_A^m

where L = k+l+m is the angular moment

In order to compute the total multipole moment, you must include the nuclear
contribution, which requires the atomic numbers and Cartesian coordinates of
the nuclei.

.. literalinclude:: ../data/examples/property/multipole_moment.py
        :lines: 1-20
        :caption: /data/examples/multipole_moment.py

Excitation energy
=================
The fundamental idea behind CIS is the representation of the excited-state wave functions as linear combinations of singly excited determinants relative to the Hartree-Fock reference wave function, viz.

.. math::
    |\Psi(m)\rangle = \sum_{jb} c_j^b(m) \{ a^+_b a_j \} |\Phi_0\rangle =
    \sum_{jb} c_j^b(m) |\phi_j^b\rangle

where m identifies the various excited states, and we will use i and j (a and b) to denote occupied (unoccupied) spin-orbitals. Inserting this into the Schrödinger equation and left-projecting onto a particular singly excited determinant gives

.. math::
    \sum_{jb} \langle \phi_i^a|H|\phi_j^b\rangle c_j^b(m) = E(m) c_i^a(m)

If we recognize that we have one of these equations for every combination of i and a spin-orbitals, then this equation may be viewed as a matrix eigenvalue problem:

.. math::
    {\mathbf {\underline{\underline H}}}\ {\mathbf {\underline c}}(m) = E(m) {\mathbf {\underline c}}(m)

.. literalinclude:: ../data/examples/property/excitation_energy.py
        :lines: 1-20
        :caption: /data/examples/excitation_energy.py
