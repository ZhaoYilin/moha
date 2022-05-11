===================
Perturbation Theory
===================
Perturbation theory is a collection of versatile methods used in many branches of science, and it divides the system into a model part :math:`\hat{H}_0` which is a known approximation to the real system :math:`\hat{H}` and a perturbation part :math:`\hat{V}`. 

.. math::

    \hat{H} = \hat{H}_0 + \lambda\hat{V}

where :math:`\lambda` is a parameter that is small enough to guarantee convergence.

Møller--Plesset perturbation theory is a particular case of perturbation theory, where we take the Fock operator :math:`\hat{F}` as the model operator, and the perturbation operator is given by

.. math::

    \hat{V} = \hat{H} - \hat{F}

As the cornerstone of ab initio quantum chemistry, Hartree--Fock calculation results in the wavefunction :math:`\Phi_{HF}`, the Fock operator :math:`\hat{F}`, ground state energy :math:`E_{HF}` and the orbital energies :math:`\{\epsilon_i\}`. Thus, the zeroth--order time--independent Schrodinger equation is defined as:
 
.. math::

    \hat{H}^{(0)}\Phi_{0}^{(0)} = E_{0}^{(0)}\Phi_{0}^{(0)}

here the zeroth--order Hamiltonian :math:`\hat{H}^{(0)}` is taken as the Fock operator，

.. math::

    \hat{H}^{(0)} = \hat{F}

and the zeroth--order wavefunction :math:`\Phi_{0}` refers to the Hartree--Fock wavefunction，

.. math::

    \Phi_{0} = \Phi_{HF}

furthermore, the corresponding energy :math:`E^{(0)}_0` is the sum of orbital energies :math:`\{\epsilon_i\}`, which is the eigenvalues of the Fock operator :math:`\hat{F}`. 

.. math::

    E^{(0)} = \sum_i^{occ} \epsilon_i

As the solution of the zeroth--order equation is known, one expects that the Møller--Plesset theory performs on top of the previous result. Therefore we can determine the various nth--order corrections.

Algorithm
=========

* **Step 1: Build The Integral**

The pre--computed quantities for the SCF algorithm include:

    * The nuclear repulsion energy :math:`E_{nuc}`.

    * One--electron integrals overlap :math:`\mathbf{S}`, kinetic energy :math:`\mathbf{T}`, and nuclear attraction :math:`\mathbf{V}`.

    * Two--electron integral, the electron repulsion :math:`(\mu\nu\vert\lambda\sigma)`.

Please check out the :ref:`Hamiltonian <Hamiltonian>` section for further information to learn more about the pre--computed quantities.

* **Step 2: Optimized the molecular orbital coefficients by SCF Calculation**

The Hartree--Fock method is an uncorrelated mean--field theory that offers a qualitative description of chemical systems. Although Hartree--Fock theory is only qualitatively correct, it forms the basis for more accurate models and becomes the cornerstone of ab initio quantum chemistry.

Please check out the :doc:`user_hf` section for further information.

* **Step 3: Transformation of atomic orbital to molecular orbital**

With the optimized LCAO_MO coefficients, we can transform the operators from the atomic orbital basis to the molecular orbital basis.

For the one electron operators:

.. math::

    h_{ij} = \sum_{\mu\nu}C_{\mu}^jC_{\nu}^i \langle\phi_{\mu}\vert\hat{h}\vert\phi_{\nu}\rangle
    = \sum_{\mu\nu} C_{\mu}^j C_{\nu}^i F_{\mu\nu}

For the two electron operators:

.. math::

    \langle pq \vert rs\rangle = \sum_{\mu} C_{\mu}^p \left[
    \sum_{\nu} C_{\nu}^q \left[
    \sum_{\lambda} C_{\lambda}^r \left[
    \sum_{\sigma} C_{\sigma}^s 
    \langle\mu\nu\vert\lambda\sigma\rangle\right]\right]\right]

* **Step 4: Calculate the Correction Energy to n--th Order**

Perturbation theory applies a biased bi--partition of the Hamiltonian :math:`\hat{H}`. The whole Hilbert space :math:`\mathcal{H}` is also split into two parts: model space :math:`\mathcal{P}` and the orthogonal space :math:`\mathcal{Q}`.

.. math::

    \mathcal{H} = \mathcal{P} \oplus \mathcal{Q}

The model space :math:`\mathcal{P}` is spanned by the reference function :math:`\{\Phi_{0}\}`, and the complementary space :math:`\mathcal{Q}` is spanned by the excitation configurations :math:`\{S, D, T, \dots\}`.

Given the Møller--Plesset perturbation operator :math:`\hat{V}`, with a reasonably accurate reference :math:`\Phi_{HF}`, we can calculate the electron correlation energy to arbitrary order in a non--iterative way. Here list the equation of correlation energy from orders one up to three:
 
* MP1

.. math::

    E_{MP1} = \langle \Phi_0 \vert\hat{V}\vert \Phi_0 \rangle
    = \frac{1}{2} \sum_{ij} \langle ij\vert\vert ij\rangle

* MP2

.. math::

    E_{MP2} = \langle \Phi_0 \vert\hat{V}\vert D \rangle
    \langle D \vert\hat{V}\vert \Phi_0 \rangle
    = \sum_{ij} \sum_{ab} 
    \frac{\langle ia\vert jb\rangle (2 \langle ia\vert jb\rangle -
    \langle ib\vert ja\rangle)}
    {\epsilon_i + \epsilon_j - \epsilon_a - \epsilon_b}

* MP3

.. math::

    E_{MP3} = \langle \Phi_0 \vert\hat{V}\vert D \rangle
    \langle D \vert\hat{V}\vert D \rangle
    \langle D \vert\hat{V}\vert \Phi_0 \rangle
    &=\frac{1}{8}\sum_{abcdrs}\frac{\left\langle ab\left|\right|rs\right\rangle \left\langle cd\left|\right|ab\right\rangle \left\langle rs\left|\right|cd\right\rangle }{\left(\epsilon_{a}+\epsilon_{b}-\epsilon_{r}-\epsilon_{s}\right)\left(\epsilon_{c}+\epsilon_{d}-\epsilon_{r}-\epsilon_{s}\right)}\\
    &+\frac{1}{8}\sum_{abrstu}\frac{\left\langle ab\left|\right|rs\right\rangle \left\langle rs\left|\right|tu\right\rangle \left\langle tu\left|\right|ab\right\rangle }{\left(\epsilon_{a}+\epsilon_{b}-\epsilon_{r}-\epsilon_{s}\right)\left(\epsilon_{a}+\epsilon_{b}-\epsilon_{t}-\epsilon_{u}\right)}\\
    &+\sum_{abcrst}\frac{\left\langle ab\left|\right|rs\right\rangle \left\langle cs\left|\right|tb\right\rangle \left\langle rt\left|\right|ac\right\rangle }{\left(\epsilon_{a}+\epsilon_{b}-\epsilon_{r}-\epsilon_{s}\right)\left(\epsilon_{a}+\epsilon_{c}-\epsilon_{r}-\epsilon_{t}\right)}

Examples
========
Perturbation theory offers robust methods for obtaining the approximate solution of problems involving small parameter :math:`\epsilon`. Moller-Plesset perturbation theory is a particular case of perturbation theory where we consider an unperturbed Hamiltonian operator :math:`\hat{H}_0`, to which we add a small perturbation :math:`\hat{V}`.

In quantum chemistry, Moller-Plesset is the most popular perturbation theory, its unperturbed Hamiltonian is the shifted Fock operator, and the zeroth-order wave function :math:`\Psi_{0}` is the lowest eigenstate of the Fock operator :math:`\hat{F}`.

* MP2

Singly excited Slater determinants do not contribute to the correction energy because of the Brillouin theorem. Hence the second-order energy is the first meaningful correction in Moller-Plesset Perturbation theory and results in the MP2 method.  

.. literalinclude:: ../../data/examples/pt/mp2.py
    :caption: /data/examples/pt/mp2.py

* MP3

After the second-order energy, the involvement of third-order correction still improves the Hartree-Fock method with cheap costs.  

.. literalinclude:: ../../data/examples/pt/mp3.py
    :caption: /data/examples/pt/mp3.py
