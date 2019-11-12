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

.. literalinclude:: ../data/examples/property/mulliken_population_analysis.py
        :lines: 1-20
        :caption: /data/examples/mulliken_population_analysis.py


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


Multipole Polarizability
========================

Second-order energy E^{(2)} is the second derivative of E(\vec{F}), is connected with the {\color{red}{polarizability}} {\color{blue}{$\alpha$}} of the molecule:

.. math::
    E^{(2)}_{ij} = -\frac{1}{2!}\alpha_{ij}

\item The \alpha_{ij} polarizability tensor is defined through the proportionality of the induced dipole moment to the field strength {\vec{F} in the limit \vec{F} = 0

.. math::
    \alpha_{ij} = \frac{\partial \mu_i}{\partial F_j}

Here we go back to the definitions of the properties as derivatives of the
energy in the presence of the perturbation by electric field.

.. math::
    \epsilon(\upsilon) = \epsilon^{(0)} +\epsilon^{(1)}\upsilon+ \frac{1}{2}\epsilon^{(2)}\upsilon^{2}+\dots
 
Rational models are defined as ratios of polynomials and are given by

.. math::
    y = \frac{\sum_{i=1}^{n+1}p_ix^{n+1-i}}{x^m+\sum_{i=1}^{m}q_ix^{m-1}}

where the coefficient associated with x^m is always 1. This
makes the numerator and denominator unique when the polynomial degrees are the same. In our polarizability calculations, the expression is:

.. math::
    E(\vec{F}) = \frac{E(0)+\sum_i b_iF_i +\sum_{ij}
    c_{ij}F_iF_j+\dots}{1+\sum_i B_iF_i} +\dots

Excitation energy
=================

Configuration Interaction Singles (CIS)
---------------------------------------
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

.. literalinclude:: ../data/examples/property/cis_excitation_energy.py
        :lines: 1-20
        :caption: /data/examples/cis_excitation_energy.py

Time-Dependent Hartree Fock (TDHF)/The Random Phase Approximation (RPA)
-----------------------------------------------------------------------
The TDHF/RPA approach is closely related to CIS in that only singles are involved in the wave function expansion, except that both excitation and “de-excitation” operators are involved. (It is probably better to view the TDHF/RPA wave function expansion in terms of orbital rotations instead of Slater determinants, but that's a discussion for another day.) The TDHF/RPA eigenvalue equations take the form

.. math::
    \left(\begin{array}{cc}
    {\mathbf {\underline{\underline A}}} & {\mathbf {\underline{\underline B}}} \\
    -{\mathbf {\underline{\underline B}}} & - {\mathbf {\underline{\underline A}}} \\
    \end{array}
    \right) \left(
    \begin{array}{c}
    {\mathbf {\underline X}} \\
    {\mathbf {\underline Y}} \\
    \end{array}
    \right) = E \left(
    \begin{array}{c}
    {\mathbf {\underline X}} \\
    {\mathbf {\underline Y}} \\
    \end{array}
    \right).

The definition of the A matrix is just the CIS matrix itself, viz.

.. math::
    A_{ia,jb} \equiv \langle \phi_i^a | H | \phi_j^b \rangle = f_{ab} \delta_{ij} - f_{ij} \delta_{ab} + \langle aj || ib \rangle,

while X and Y are the parameters of single excitations and de-excitations, respectively, and the B matrix is simply

.. math::
    B_{ia,jb} \equiv \langle ab || ij \rangle.

.. literalinclude:: ../data/examples/property/tdhf_excitation_energy.py
        :lines: 1-20
        :caption: /data/examples/tdhf_excitation_energy.py
