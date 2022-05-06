===================
Perturbation Theory
===================
Perturbation theory offers robust methods for obtaining the approximate solution of problems involving small parameter :math:`\epsilon`. Moller-Plesset perturbation theory is a particular case of perturbation theory where we consider an unperturbed Hamiltonian operator :math:`\hat{H}_0`, to which we add a small perturbation :math:`\hat{V}`.

In quantum chemistry, Moller-Plesset is the most popular perturbation theory, its unperturbed Hamiltonian is the shifted Fock operator, and the zeroth-order wave function :math:`\Psi_{0}` is the lowest eigenstate of the Fock operator :math:`\hat{F}`.


For perturbation theory, we divide the Hamiltonian into zero--order operator :math:`\hat{H}_0` and perturbaing operator :math:`\lambda\hat{V}`. Specifically in M{\o}ller--Plesset theory, we use the Fock operator as the zero--order operator and the Hartree--Fock state as the zero--order state:

.. math::

    \hat{H}_0 = \hat{F}; 
    \quad \hat{V} = \hat{H} - \hat{F} 

the unperturbed Schr\"{o}dinger equation:

.. math::

    \quad \hat{F}|\Psi_{HF}\rangle = \sum_{i}\varepsilon_i \vert\Psi_{HF}\rangle; 

the time-independent related Schr\"{o}dinger equation:

.. math::

    (\hat{F} + \lambda \hat{V})|\Psi_n \rangle = E_n \vert\Psi_n \rangle

we apply Taylor series to expands the energies and states in eqn \ref{eq:mp3} of the parameter $\lambda$

.. math::

    E_n = E_n^{(0)} + \lambda E_n^{(1)} + \lambda^2 E_n^{(2)} + \dots

and

.. math::

    \vert\Psi_{n}\rangle = \vert\Psi_{n}^{(0)} \rangle 
    + \lambda \vert\Psi_{n}^{(0)} \rangle + \lambda^2 \vert\Psi_{n}^{(2)} \rangle + \dots

the :math:`E_n^{(k)}` and :math:`\Psi_n^{(k)}` are the corresponding derivatives

.. math::

    E_n^{(k)} = \frac{1}{k!} \frac{\partial^k E_n}{\partial \lambda^k}

.. math::

    \Psi_n^{(k)} = \frac{1}{k!} \frac{\partial^k \Psi_n}{\partial \lambda^k}

substitute of these expand series into the time-independent Schr\"{o}dinger equation

.. math::

    (\hat{F} + \lambda \hat{V})\vert\Psi_n^{(0)} + \lambda \vert\Psi_n^{(1)} + \lambda^2 \vert\Psi_n^{(2)} + \dots\rangle&= \\
    (E_n^{(0)} + \lambda E_n^{(1)} + \lambda^2 E_n^{(2)} +\dots) \vert\Psi_n^{(0)} + \lambda \vert\Psi_n^{(1)} + \lambda^2 \vert\Psi_n^{(2)} + \dots\rangle

we expand both side of the equation, put the terms involving equate powers of :math:`\lambda` and generate a set
of equations. 

.. math::

    \hat{F}|\Psi_n^{(0)}\rangle 
    &= E_n^{(0)}|\Psi_n^{(0)}\rangle \\
    \hat{F}|\Psi_n^{1}\rangle +\hat{V}|\Psi_n^{(0)}\rangle	
    &= E_n^{(0)}|\Psi_n^{(1)}\rangle + E_n^{(1)}|\Psi_n^{(0)}\rangle\\
    \hat{F}|\Psi_n^{2}\rangle +\hat{V}|\Psi_n^{(1)}\rangle 	
    &= E_n^{(0)}|\Psi_n^{(2)}\rangle + E_n^{(1)}|\Psi_n^{(1)}\rangle + E_n^{(2)}|\Psi_n^{(0)}\rangle\\
    &\dots				

as the solution of zeroth--order equation is known, one expects, Möller--Plesset theory approaches are performed on top of previous result. Therefore we can determine the various nth--order corrections.

Given a reasonably accurate Hartree--Fock reference, we may improve on it by Möller--Plesset perturbation theory. As a non-iterative way to calculate electron correlation energy, The MP2 energy represents a highly successful approximation at a fraction of the cost compare with other post--Hartree--Fock methods. 

Unfortunately in a sufficiently large system, the M{\o}ller--Plesset method often diverges. Higher-order corrections may be calculated, but convergence is often poor. 

Algorithm
=========

* **Step 1: Build The Integral**

The pre--computed quantities for the SCF algorithm include:

    * The nuclear repulsion energy :math:`E_{nuc}`.

    * One--electron integrals overlap :math:`\mathbf{S}`, kinetic energy :math:`\mathbf{T}`, and nuclear attraction :math:`\mathbf{V}`.

    * Two--electron integral, the electron repulsion :math:`(\mu\nu\vert\lambda\sigma)`.

Please check out the :ref:`Hamiltonian <Hamiltonian>` section for further information to learn more about the pre--computed quantities.

* **Step 2: Optimized the molecular orbital coefficients by SCF Calculation**

The Hartree-Fock method is an uncorrelated mean-field theory that offers a qualitative description of chemical systems. Although Hartree-Fock theory is only qualitatively correct, it forms the basis for more accurate models and becomes the cornerstone of ab initio quantum chemistry.

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

MP2
---
Moller-Plesset Second Order Perturbation theory (MP2) is a non-iterative way to calculate electron correlation energy. The expression for it's energy contribution can be written as:

.. math::

    E_{\rm MP2} = \sum_{ij} \sum_{ab} 
    \frac{\langle ia\vert jb\rangle (2 \langle ia\vert jb\rangle -
    \langle ib\vert ja\rangle)}
    {\epsilon_i + \epsilon_j - \epsilon_a - \epsilon_b}

To run a MP2 calculation in MoHa, the example is below:

.. literalinclude:: ../../data/examples/pt/mp2.py
    :caption: /data/examples/pt/mp2.py

MP3
---
The third order Møller-Plesset correction to the energy is given as:

.. math::

    E_{MP3}&=\frac{1}{8}\sum_{abcdrs}\frac{\left\langle ab\left|\right|rs\right\rangle \left\langle cd\left|\right|ab\right\rangle \left\langle rs\left|\right|cd\right\rangle }{\left(\epsilon_{a}+\epsilon_{b}-\epsilon_{r}-\epsilon_{s}\right)\left(\epsilon_{c}+\epsilon_{d}-\epsilon_{r}-\epsilon_{s}\right)}\\
    &+\frac{1}{8}\sum_{abrstu}\frac{\left\langle ab\left|\right|rs\right\rangle \left\langle rs\left|\right|tu\right\rangle \left\langle tu\left|\right|ab\right\rangle }{\left(\epsilon_{a}+\epsilon_{b}-\epsilon_{r}-\epsilon_{s}\right)\left(\epsilon_{a}+\epsilon_{b}-\epsilon_{t}-\epsilon_{u}\right)}\\
    &+\sum_{abcrst}\frac{\left\langle ab\left|\right|rs\right\rangle \left\langle cs\left|\right|tb\right\rangle \left\langle rt\left|\right|ac\right\rangle }{\left(\epsilon_{a}+\epsilon_{b}-\epsilon_{r}-\epsilon_{s}\right)\left(\epsilon_{a}+\epsilon_{c}-\epsilon_{r}-\epsilon_{t}\right)}

To run a MP2 calculation in MoHa, the example is below:

.. literalinclude:: ../../data/examples/pt/mp3.py
    :caption: /data/examples/pt/mp3.py
