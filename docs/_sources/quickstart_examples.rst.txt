========
Examples
========
This chapter introduces the primary use of MoHa calculation, and you can find those examples in the :file:`moha/data/examples` directory. The accompanying user guide will present additional details on the various modules.

Hartree-Fock
============
The Hartree-Fock method is an uncorrelated mean-field theory that offers a qualitative description of chemical systems. Although Hartree-Fock theory is only qualitatively correct, it forms the basis for more accurate models and becomes the cornerstone of ab initio quantum chemistry.

Plain Solver
------------
A basic Hartree-Fock calculation with a plain solver can be performed by:

.. literalinclude:: ../../data/examples/hf/scf.py
    :caption: /data/examples/hf/scf.py

DIIS Solver
------------
The convergence of the plain Hartree-Fock solver is relatively slow. One of the most effective methods for speeding up such iterative sequences is the "Direct Inversion in the Iterative Subspace(DIIS)" approach.

.. literalinclude:: ../../data/examples/hf/scf_diis.py
    :caption: /data/examples/hf/scf_diis.py

Post-Hartree-Fock
=================
Quantum chemistry is mature enough to have canonical approaches followed by the Hartree-Fock method. Their behavior, successes, and most importantly, their failures are widely known and accepted. With a couple of exceptions which we will not mention, they are:

- Configuration Interaction(CI)
- Coupled-Cluster(CC)
- Perturbation Theory(PT)

Configuration Interaction
-------------------------
Configuration interaction approximates wave function by a linear expansion of N-electron basis functions made of a given one-electron basis. With CI wave function and CI basis set, the Schrödinger equation becomes a matrix-eigenvalue equation.

* Full CI

The full configuration interaction(FCI) method assumes that all electrons are correlated among all orbitals in a given system. Hence it provides numerically exact solutions (within the infinitely flexible complete basis set) to the electronic time-independent, non-relativistic Schrödinger equation.

.. literalinclude:: ../../data/examples/ci/fci.py
    :caption: /data/examples/ci/fci.py

* CISD

FCI method is exact in a given atomic orbital basis but prohibitively expensive. Therefore, to balance accuracy and computational time, we truncate the CI space by excitation level relative to the reference state. The most common truncation of the CI space expansion is CI singles and doubles (CISD).

.. literalinclude:: ../../data/examples/ci/cisd.py
    :caption: /data/examples/ci/cisd.py

Coupled-Cluster
---------------
Coupled cluster methods are among the most accurate electronic structure methods available today. Its wave function ansatz gives in exponential form :math:`e^{\hat{T}} \vert\Phi_0\rangle`

Where the cluster operator :math:`\hat{T}` is the sum of :math:`\{\hat{T}_{i}\}` generates all possible determinants having the i-th excitations from the reference :math:`\vert\Phi_0\rangle`.

* CCSD

Truncating the cluster operator :math:`\hat{T}` at doubles leads to the coupled cluster singles and doubles method(CCSD). By keeping only :math:`\hat{T}_1` and :math:`\hat{T}_2`, we reach the balance between accuracy and computational time.

.. literalinclude:: ../../data/examples/cc/ccsd.py
    :caption: /data/examples/cc/ccsd.py

* CCSD(T)

Based on CCSD, the addition of perturbative triple correction gives the CCSD(T) method the “gold standard” in computational chemistry. It is now one of the most accurate methods applicable to reasonably large molecules.

.. literalinclude:: ../../data/examples/cc/ccsd_t.py
    :caption: /data/examples/cc/ccsd_t.py

Perturbation Theory
-------------------
Perturbation theory is a collection of versatile methods used in many branches of science. It divides the system into a model part :math:`\hat{H}_0` which is a known approximation to the real system :math:`\hat{H}` and a perturbation part :math:`\lambda\hat{V}`, where :math:`\lambda` is a parameter that is small enough to guarantee convergence.

Møller--Plesset perturbation theory is a particular case of perturbation theory, where we take the Fock operator :math:`\hat{F}` as the model operator.

* MP2

Depending on the choice of Møller--Plesset perturbation operator :math:`\hat{V}`, the first order correction energy :math:`E_{MP1}` might be zero or none zero, but :math:`E_{MP0} + E_{MP1}` is always equal to the Hartree--Fock energy :math:`E_{HF}`. Hence, the first meaningful correction in Møller--Plesset perturbation theory is second-order energy.

.. literalinclude:: ../../data/examples/pt/mp2.py
    :caption: /data/examples/pt/mp2.py

* MP3

After the second-order energy, the involvement of third-order correction still improves the Hartree-Fock method with cheap costs.  

.. literalinclude:: ../../data/examples/pt/mp3.py
    :caption: /data/examples/pt/mp3.py


Property
========
This section contains information about atomic and molecular properties implemented in MoHa.

Population Analyasis
--------------------
Population analysis is the study of charge distribution in molecules. It requires the overlap integrals and the electron density and information about the number of basis functions centered on each atom.

.. literalinclude:: ../../data/examples/property/population_analysis.py
    :caption: /data/examples/population_analysis.py

Multipole Moment
----------------
Multipole moments are the coefficients of a series expansion of a potential due to continuous or discrete sources. A multipole moment usually involves powers (or inverse powers) of the distance to the origin and some angular dependence.

.. literalinclude:: ../../data/examples/property/multipole_moment.py
    :caption: /data/examples/multipole_moment.py

Excitation Energy
-----------------
The simplest ab initio method to calculate excited electronic states is configuration interaction singles (CIS). The fundamental idea behind CIS is the representation of the excited-state wave functions as linear combinations of singly excited determinants relative to the Hartree-Fock reference wave function.

.. literalinclude:: ../../data/examples/property/excitation_energy.py
    :caption: /data/examples/excitation_energy.py
