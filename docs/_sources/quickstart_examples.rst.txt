========
Examples
========
This chapter provides an introduction to the basic use of MoHa calculation, you can find those examples in the :file:`moha/data/examples` directory. Additional details on the various modules are presented in the accompanying user guide.

Hartree-Fock
============
The Hartree-Fock method provides the fundamental approximation of wave function theory. Although Hartree-Fock theory is only qualitatively correct, it forms the basis for more accurate models and becomes the cornerstone of ab initio quantum chemistry.

The Hartree-Fock method assumes that the ground state wave function is the direct product state of a series of single electron wave functions. The Hartree-Fock wave function ansatz can be written as

.. math:: \vert {\kappa}\rangle = exp(-\hat{\kappa})\vert 0 \rangle

Plain Solver
------------
A basic Hartree-Fock calculation with a plain solver can be performed by:

.. literalinclude:: ../../data/examples/hf/scf.py
    :caption: /data/examples/hf/scf.py

DIIS Solver
------------
The convergence of the plain Hartree-Fock solver is rather slow, One of the most effective methods for speeding up such iterative sequences is the "Direct Inversion in the Iterative Subspace(DIIS)" approach, which was presented by Pulay in the 1980s.

.. literalinclude:: ../../data/examples/hf/scf_diis.py
    :caption: /data/examples/hf/scf_diis.py

Post-Hartree-Fock
=================
Quantum chemistry is mature enough to have canonical approaches followed by Hartree-Fock method. Their behavior, successes, and most importantly their failures are widely known and accepted. With a couple of exceptions which we will not mention, they are known as:

- Configuration Interaction(CI)
- Coupled-Cluster(CC)
- Perturbation Theory(PT)

Configuration Interaction
-------------------------
The configuration interaction method is a commonly used quantum chemistry method, which is used to correct electron correlation based on the Hartree-Fock method. This method writes the wave function as a linear combination of configuration state functions.

.. math:: \vert \Psi_{CI} \rangle = \sum_{\lbrace i \rbrace} c_i \vert i \rangle

* Full CI

The full CI method assumes that all electrons of a given system are correlated among all orbitals. Hence it provides numerically exact solutions (within the infinitely flexible complete basis set) to the electronic time-independent, non-relativistic Schrödinger equation.  

.. math::

    \vert \Psi_{FCI} \rangle = c_0 \vert \Phi_0\rangle
    + \mathop{\sum_{\lbrace i \rbrace}}_{\lbrace a \rbrace} c_i^a \vert\Phi_i^a\rangle
    + \mathop{\sum_{\lbrace i,j \rbrace}}_{\lbrace a,b \rbrace} c_{ij}^{ab} \vert\Phi_{ij}^{ab}\rangle
    + \dots

.. literalinclude:: ../../data/examples/ci/fci.py
    :caption: /data/examples/ci/fci.py


* CISD

In practice, we truncated the CI space by excitation order to save computational time. Configuration interaction with singly and doubly excited determinants methods(CISD) are implemented in many standard packages.

.. math::

    \vert\Psi_{CISD}\rangle = c_0 \vert\Phi_0\rangle
    + \mathop{\sum_{\lbrace i \rbrace}}_{\lbrace a \rbrace} c_i^a \vert\Phi_i^a\rangle
    + \mathop{\sum_{\lbrace i,j \rbrace}}_{\lbrace a,b \rbrace} c_{ij}^{ab} \vert\Phi_{ij}^{ab}\rangle

.. literalinclude:: ../../data/examples/ci/cisd.py
    :caption: /data/examples/ci/cisd.py

Coupled-Cluster
---------------
Coupled cluster methods are among the most accurate electronic structure methods available today. Its wave function ansatz is usually written in exponential form:

.. math:: \vert\Psi\rangle = e^{\hat{T}} \vert\Phi_0\rangle

Where the cluster operator :math:`\hat{T}` is defined as the sum of :math:`\hat{T}_{i}` generates all possible determinants having the i-th excitations from the reference.

.. math:: \hat{T} = \sum_{1}^{N} \hat{T}_{i}

* CCSD

In practice, we only keep the first two cluster operators :math:`\hat{T}_1` and :math:`\hat{T}_2` to reach the balance between accuracy and computational time.

To run a CCSD calculation in MoHa, the example is below:

.. literalinclude:: ../../data/examples/cc/ccsd.py
    :caption: /data/examples/cc/ccsd.py

* CCSD(T)

The CCSD(T) method is often called the "gold standard" of computational chemistry because it is one of the most accurate methods applicable to reasonably large molecules. To find the perturbative triples correction, the disconnected and connected T3 have to be calculated.

To run a CCSD(T) calculation, the example is below:

.. literalinclude:: ../../data/examples/cc/ccsd_t.py
    :caption: /data/examples/cc/ccsd_t.py

Perturbation Theory
-------------------
Perturbation theory offers a set of powerful methods for obtaining the approximate solution of problems involving small parameter :math:`\epsilon`. Moller-Plesset perturbation theory is a special case of perturbation theory, where we consider an unperturbed Hamiltonian operator :math:`\hat{H}_0`, to which a small perturbation :math:`\hat{V}` is added:

.. math::

    \hat{H} = \hat{H}_0 + \lambda\hat{V}

Here, λ is an arbitrary real parameter that controls the size of the perturbation. 

In quantum chemistry, Moller-Plesset is the most popular perturbation theory, its unperturbed Hamiltonian is defined as the shifted Fock operator, 

.. math::

   \hat{H}_0 = \hat{F} + \langle \Psi_{0} \vert (\hat{H} - \hat{F})\vert \Psi_{0} \rangle 

where the zeroth-order wave function :math:`\Psi_{0}` is the lowest eigenstate of the Fock operator :math:`\hat{F}`.

* MP2

Singly excited Slater determinants do not contribute to the correction energy because of the Brillouin theorem, the second-order energy is the first meaningful correction that appears in Moller-Plesset Perturbation theory.

To run a Moller-Plesset perturbation theory(MP2) calculation in MoHa, the example is below:

.. literalinclude:: ../../data/examples/pt/mp2.py
    :caption: /data/examples/pt/mp2.py

* MP3

After the second-order energy, the involvement of third-order correction still offers an improvement on the Hartree-Fock method with cheap costs.  

To run a MP3 calculation in MoHa, the example is below:
   
.. literalinclude:: ../../data/examples/pt/mp3.py
    :caption: /data/examples/pt/mp3.py


Property
========
This section containts information about atomic and molecular properties that can be calculated.

Population Analyasis
--------------------
Population analysis is the study of charge distribution in molecules, it requires the overlap integrals and the electron density, in addition to information about the number of basis functions centered on each atom. The charge on atom A may be computed as:

.. math:: q_A = Z_A - 2 \sum_{\mu \in A} ({\mathbf D S})_{\mu\mu},

where the summation is limited to only those basis functions centered on atom A.

.. literalinclude:: ../../data/examples/property/population_analysis.py
    :caption: /data/examples/population_analysis.py

Multipole Moment
----------------
The calculation of one-electron properties requires density a matrix and the relevant property integrals. The electronic contribution to the electric multipole moment may be computed using:

.. math:: {\mu}^{klm}=-\sum_{i}\langle \Psi_0|x_i^ky_i^lz_i^m|\Psi_0 \rangle+\sum_A Z_Ax_A^ky_A^lz_A^m

where L = k+l+m is the angular moment

To compute the total multipole moment, you must include the nuclear contribution, which requires the atomic numbers and Cartesian coordinates of the nuclei.

.. literalinclude:: ../../data/examples/property/multipole_moment.py
    :caption: /data/examples/multipole_moment.py

Excitation Energy
-----------------  
The excitation energies can be calculated by using the Configuration Interaction Singles(CIS) method. The fundamental idea behind CIS is the representation of the excited-state wave functions as linear combinations of singly excited determinants relative to the Hartree-Fock reference wave function.

.. math::

    \vert\Psi_{CIS}\rangle = c_0 \vert\Phi_0\rangle
    + \mathop{\sum_{\lbrace i \rbrace}}_{\lbrace a \rbrace} c_i^a \vert\Phi_i^a\rangle

To run an excitation energies calculation in MoHa, the example is below:

.. literalinclude:: ../../data/examples/property/excitation_energy.py
    :caption: /data/examples/excitation_energy.py
