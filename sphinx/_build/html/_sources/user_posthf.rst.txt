Post Hartree-Fock
#################

Quantum chemistry is mature enough to have a canonical approaches followed by
Hartree-Fock method. Their behaviour, successes, and most importantly their
failures are widely known and accepted. With a couple of exceptions which we
will not mention, they are known as:

- Configuration Interaction(CI)
- Coupled-Cluster(CC)
- Perturbation Theory(PT)

Configuration Interaction
=========================
CID
---
TODO

CISD
----
TODO

CISD(T)
-------
TODO

Full CI
-------
TODO

Coupled-Cluster
===============
Coupled cluster methods are among the most accurate electronic structure
methods available today. For example, with a good choice for a basis, the
CCSD(T) equations will give you the correct energy of a molecular system to
chemical accuracy (~0.1 kcal/mol). 

CCD
---
TODO

CCSD
----
To implement the  CCSD method, we should follow the procedure:

* Preparing the Spin-Orbital Basis Integrals
* Build the Initial-Guess Cluster Amplitudes
For Hartree-Fock reference determinants, the most common initial guess for the
cluster amplitudes are the Moller-Plesset first-order perturbed wave function:

.. math::
    t_i^a = 0
.. math::
    t_{ij}^{ab} = \frac{\langle ij || ab \rangle}{\epsilon_i + \epsilon_j -
    \epsilon_a - \epsilon_b}

* Calculate the CC Intermediates
Use the spin-orbital Eqs. 3-13 from Stanton's paper to build the two-index (F)
and four-index (W) intermediates, as well as the effective doubles (labelled
with the Greek letter tau). 

* Compute the Updated Cluster Amplitudes
Use Eqs. 1 and 2 from Stanton's paper to compute the updated T1 and T2 cluster amplitudes.

* Check for Convergence and Iterate Calculate the current CC correlation energy:

.. math::
    E_{\rm CC} = \sum_{ia} f_{ia} t_i^a + \frac{1}{4} \sum_{ijab} \langle
    ij||ab\rangle t_{ij}^{ab} + \frac{1}{2} \sum_{ijab} \langle ij||ab\rangle t_i^a
    t_j^b

Compare energies and cluster amplitudes between
iterations to check for convergence to some specified cutoff. If convergence is
reached, you're done; if not, return to third step and continue.

To run a CCSD calculation in MoHa, the example is below:

.. literalinclude:: ../data/examples/posthf/ccsd.py
        :lines: 1-20
        :caption: /data/examples/posthf/ccsd.py

Perturbation Theory
===================

MP2
---
Moller-Plesset Second Order Perturbation theory (MP2) is a non-iterative
way to calculate electron correlation energy The expression for it's energy
contribution can be written as:

.. math::
    E_{\rm MP2} = \sum_{ij} \sum_{ab} \frac{(ia|jb) \left[ 2 (ia|jb) -
    (ib|ja) \right]}{\epsilon_i + \epsilon_j - \epsilon_a - \epsilon_b}

where the two electron integral are molecular integral and should be transform
from atomic orbtial to molecular orbital basis.

.. math::
    (pq|rs) = \sum_\mu C_\mu^p \left[ \sum_\nu C_\nu^q \left[ \sum_\lambda
    C_\lambda^r \left[ \sum_\sigma C_\sigma^s (\mu \nu|\lambda \sigma) \right]
    \right] \right]

To run a MP2 calculation in MoHa, the example is below:

.. literalinclude:: ../data/examples/posthf/mp2.py
        :lines: 1-20
        :caption: /data/examples/posthf/mp2.py

MP3
---
TODO 
