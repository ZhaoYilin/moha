===============
Coupled-Cluster
===============
Coupled cluster methods are among the most accurate electronic structure methods available today. Its wave function ansatz is usually written in exponential form:

.. math:: \vert\Psi\rangle = e^{\hat{T}} \vert\Phi_0\rangle

Where the cluster operator :math:`\hat{T}` is defined as the sum of :math:`\hat{T}_{i}` generates all possible determinants having the i-th excitations from the reference.

.. math:: \hat{T} = \sum_{i=1}^{N} \hat{T}_{i}

the first few components of cluster operator are

.. math::

    \hat{T}_1 &= \mathop{\sum_{i\in occ}}_{\alpha \in virt} t_i^{\alpha} a_{\alpha}^{\dagger} a_{i}\\
    \hat{T}_2 &= \frac{1}{(2!)^2}\mathop{\sum_{ij\in occ}}_{\alpha \beta \in virt} 
    t_{ij}^{\alpha\beta}a_{\alpha}^{\dagger} a_{\beta}^{\dagger} a_{j}a_{i}

apply a similarity transformation, multiplying the CC Schrödinger equation in the form

.. math::

    \hat{H} e^{\hat{T}} \vert\Phi_{0} \rangle &= E e^{\hat{T}} \vert\Phi_{0} \rangle \\
    e^{\hat{-T}}\hat{H} e^{\hat{T}} \vert\Phi_{0} \rangle &= e^{\hat{-T}}E e^{\hat{T}} \vert\Phi_{0} \rangle \\
    e^{\hat{-T}}\hat{H} e^{\hat{T}} \vert\Phi_{0} \rangle &= E \vert\Phi_{0} \rangle \\

the amplitudes should be optimized to determine the CC wave function. Unlike Hartree--Fock or CI wave function is optimized by variational principle, CC wave function is optimized by projecting. 

.. math::

    \langle \Phi_0 \vert\hat{H} e^{\hat{T}} \vert\Phi_{0} \rangle &= E\\
    \langle \Phi_{ij\dots}^{ab\dots}\vert e^{\hat{-T}}\hat{H} e^{\hat{T}} \vert\Psi_{0} \rangle &= 0 \\

by solving the nonlinear equations, the CC energy and amplitudes are determined.
 
Algorithm
=========

* **Step 1: Build the Integral**

The pre--computed quantities for the SCF algorithm include:

    * The nuclear repulsion energy :math:`E_{nuc}`.

    * One--electron integrals overlap :math:`\mathbf{S}`, kinetic energy :math:`\mathbf{T}`, and nuclear attraction :math:`\mathbf{V}`.

    * Two--electron integral, the electron repulsion :math:`(\mu\nu\vert\lambda\sigma)`.

Please check out the :ref:`Hamiltonian <Hamiltonian>` section for further information to learn more about the pre--computed quantities.

* **Step 2: Optimization of Molecular Orbital Coefficients**

The Hartree-Fock method is an uncorrelated mean-field theory that offers a qualitative description of chemical systems. Although Hartree-Fock theory is only qualitatively correct, it forms the basis for more accurate models and becomes the cornerstone of ab initio quantum chemistry.

Please check out the :doc:`user_hf` section for further information.

* **Step 3: Transformation of Atomic Orbital to Molecular Orbital**

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

* **Step 4: Build the Initial-Guess Cluster Amplitudes**

For Hartree-Fock reference determinants, the most common initial guess for the cluster amplitudes are the Moller-Plesset first-order perturbed wave function:

.. math::
    t_i^a = 0

.. math::
    t_{ij}^{ab} = \frac{\langle ij \vert\vert ab \rangle}{\epsilon_i + \epsilon_j - \epsilon_a - \epsilon_b}

* **Step 5: Calculate the CC Intermediates**

The definition of three two--index intermediates :math:`\mathcal{F}`:

.. math::

    \mathcal{F}_{ae} &= (1-\delta_{ae})f_{ae} -\frac{1}{2} \sum_m f_{me}t_m^a\\
    &+\sum_{mf}t_m^f \langle ma\vert\vert fe\rangle
    -\frac{1}{2}\sum_{mnf}\tilde{\tau}_{mn}^{af} \langle mn\vert\vert ef\rangle


.. math::

    \mathcal{F}_{mi} &= (1-\delta_{mi})f_{mi} -\frac{1}{2} \sum_e t_i^e f_{me}\\
    &+\sum_{en}t_n^e \langle mn \vert\vert ie\rangle
    +\frac{1}{2}\sum_{nef}\tilde{\tau}_{in}^{ef} \langle mn\vert\vert ef\rangle


.. math::

    \mathcal{F}_{me} = f_{me} + \sum_{nf} t_n^f \langle mn\vert\vert ef\rangle

The definition of three four--index intermediates :math:`\mathcal{W}`:

.. math::

    \mathcal{W}_{mnij} &= \langle mn \vert\vert ij \rangle 
    + P_{(ij)} \sum_e t_j^e \langle mn\vert\vert ie \rangle\\
    &+\frac{1}{4} \sum_{ef}\tau_{ij}^{ef} \langle mn \vert\vert ef\rangle

.. math::

    \mathcal{W}_{abef} &= \langle ab \vert\vert ef \rangle 
    - P_{(ij)} \sum_m t_m^b \langle am\vert\vert ef \rangle\\
    &+\frac{1}{4} \sum_{mn}\tau_{mn}^{ab} \langle mn \vert\vert ab\rangle

.. math::

    \mathcal{W}_{mbej} &= \langle mb \vert\vert ej \rangle 
    +\sum_f t_j^f \langle mb \vert\vert ef \rangle\\
    &-\sum_{n}t_n^b \langle mn \vert\vert ej\rangle
    -\sum_{nf} (\frac{1}{2}t_{jn}^{fb}+t_j^ft_n^b)\langle mn\vert\vert ef\rangle

The definition of effective double excitation operaotrs :math:`\tau` and :math:`\tilde{\tau}` 

.. math::

    \tilde{\tau}_{ij}^{ab} &= t_{ij}^{ab} +\frac{1}{2} (t_i^at_j^b-t_i^bt_j^a)\\
    \tau_{ij}^{ab} &= t_{ij}^{ab} + t_i^at_j^b-t_i^bt_j^a

The denominator arrays are defined by

.. math::

    D_i^a &= f_{ii} - f_{aa}\\
    D_{ij}^{ab} &= f_{ii} + f_{jj} - f_{aa} - f_{bb}\\




* **Step 6: Compute the Updated Cluster Amplitudes**

Equation to compute the T1 cluster amplitudes.

.. math::

    t_i^aD_i^a &= f_{ia} + \sum_e t_i^e \mathcal{F}_{ae} 
    - \sum_m t_m^a\mathcal{F}_mi +\sum_{me}t_{im}^{ae}\mathcal{me}\\
    &- \sum_{nf} t_n^f \langle na\vert\vert if \rangle 
    -\frac{1}{2}\sum_{mef}t_{im}^{ef}\langle ma\vert\vert ef\rangle\\
    &-\frac{1}{2}\sum_{men}t_{mn}^{ae}\langle nm\vert\vert ei\rangle\\

Equation to compute the T2 cluster amplitudes.

.. math::

    t_{ij}^{ab}D_{ij}^{ab} &= \langle ij\vert\vert ab\rangle 
    +P_{(ab)}\sum_e t_{ij}^{ae}(\mathcal{F}_{be}-\frac{1}{2} \sum_m t_m^b\mathcal{F}_{me})\\
    &-P_{(ij)}\sum_m t_{im}^{ab}(\mathcal{F}_{mj}+\frac{1}{2} \sum_e t_j^e\mathcal{F}_{me})\\
    &+\frac{1}{2}\sum_{mn}\tau_{mn}^{ab}\mathcal{W}_{mnij} 
    +\frac{1}{2} \sum_{ef}\tau_{ij}^{ef}\mathcal{W}_{abef}\\
    &+P_{(ij)}P_{(ab)}\sum_{me}(t_{im}^{ae}\mathcal{W}_{mbej}
    -t_i^et_m^a\langle mb\vert\vert ej\rangle)\\
    &+P_{(ij)}\sum_{e}t_{i}^{e}\langle ab\vert\vert ej\rangle
    -P_{(ab)} \sum_{m}t_m^a \langle mb\vert\vert ij\rangle\\


Check for Convergence and Iterate Calculate the current CC correlation energy:

.. math::
    E_{\rm CC} = \sum_{ia} f_{ia} t_i^a + \frac{1}{4} \sum_{ijab} \langle
    ij||ab\rangle t_{ij}^{ab} + \frac{1}{2} \sum_{ijab} \langle ij||ab\rangle t_i^a
    t_j^b

Compare energies and cluster amplitudes between iterations to check for convergence to some specified cutoff. If convergence is reached, you're done; if not, return to third step and continue.

* **Step 7: Find the Perturbative Triples Correction**

the disconnected and connected T3 have to be calculated. The disconnected is found as:

.. math::

    D^{abc}_{ijk}t^{abc}_{ijk,disconnected}=P(i/jk)P(a/bc)t^{a}_{i}\langle jk \vert\vert bc\rangle

And the connected is found as:

.. math::

    D_{ijk}^{abc}t_{ijk,\mathrm{connected}}^{abc}=P\left(i/jk\right)P\left(a/bc\right)\left[\sum_{e}^{virt}t_{jk}^{ae}\left\langle ei\left|\right|bc\right\rangle -\sum_{m}^{occ}t_{im}^{bc}\left\langle ma\left|\right|jk\right\rangle \right]

In the above equations the following definitions is used:

.. math::

    D_{ijk}^{abc}=f_{ii}+f_{jj}+f_{kk}-f_{aa}-f_{bb}-f_{cc}

.. math::

    P\left(i/jk\right)f\left(i,j,k\right)=f\left(i,j,k\right)-f\left(j,i,k\right)-f\left(k,j,i\right)

The energy correction can now be found as:

.. math::

    E_{\mathrm{\left(T\right)}}=\frac{1}{36}\sum_{i,j,k,a,b,c}^{occ,occ,occ,virt,virt,virt}t_{ijk,\mathrm{connected}}^{abc}D_{ijk}^{abc}\left(t_{ijk,\mathrm{connected}}^{abc}+t_{ijk,\mathrm{disconnected}}^{abc}\right)

Examples
========   
* CCSD

Truncating the cluster operator :math:`\hat{T}` at doubles leads to the coupled cluster singles and doubles method(CCSD). By keeping only :math:`\hat{T}_1` and :math:`\hat{T}_2`, we reach the balance between accuracy and computational time.

.. literalinclude:: ../../data/examples/cc/ccsd.py
    :caption: /data/examples/cc/ccsd.py

* CCSD(T)

Based on CCSD, the addition of perturbative triple correction gives the CCSD(T) method the “gold standard” in computational chemistry. It is now one of the most accurate methods applicable to reasonably large molecules.

.. literalinclude:: ../../data/examples/cc/ccsd_t.py
    :caption: /data/examples/cc/ccsd_t.py
