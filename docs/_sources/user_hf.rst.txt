Hartree-Fock
############
Hartree Fock method provides the fundamental approximation of wave--function theory. Although Hartree Fock theory is only qualitatively correct, it forms the basis for more accurate models and become the cornerstone of ab initio quantum chemistry.

Hartree--Fock method assumes that the ground state wave function is the direct product state of a series of single electron wave--functions. The Hartree--Fock wave--function ansatz can be written as

.. math::

    \vert \kappa \rangle = exp(-\hat{\kappa})\vert 0\rangle

where the anti--Hermitian operator :math:`\hat{\kappa}` in general form is 

.. math::

    \hat{\kappa} = \sum_{pq} \kappa_{PQ} a_P^{\dagger}a_Q

The Hartree--Fock wave--function is optimized with respect of variations of spin orbitals.

.. math::

    E_{HF} = min \langle {\kappa} \vert \hat{H} \vert {\kappa} \rangle

a new set of spin orbitals can be update by unitary transformation

.. math::

    \tilde{\phi}_P = \sum_{Q} \phi_Q U_{QP}

the unitary matrix can be written as

.. math::

    {U} = exp(-{\kappa})

the spin orbitals is determined by solving a set of effective one--electron Schr\"{o}dinger equations, where the effective Hamiltonian associated is Fock operator. The eigenvalues of the Fock operator are the spin energies and eigenvectors the spin orbitals.

.. math::

    \hat{f} = \hat{h} + \hat{V}

:math:`\hat{h}` is the one--body term in chemical Hamiltonian, the two--body potential terms in original Hamiltonian is replaced by the effective one electron Fock potential :math:`\hat{V}`.

.. math::

    \hat{V} = \sum_{pq} \sum_{i} (2g_{pqii}-g_{piiq})E_{pq}

since Fock matrix is define by its own eigenvectors, an iterative procedure, self--consistent--field(SCF) method is applied. Optimized by variational theorem, Hartree--Fock energy is an upper bound to the exact energy. It typical gives errors of :math:`0.5\%` in the energy, :math:`1\%`  in molecular geometry and :math:`5-10\%` in other properties. The difference between the Hartree--Fock limit energy and the exact non--relativistic energy of a system is defined as correlation energy.

.. math::

    E_{corr} = \mathcal{E}_0 - E_{HF}


Plain Solver
============
For the Roothan-Hartee-Fock equations an orthogonal basis is needed, first the first orthogonalization matrix is constructed from the overlap matrix. First by a diagonalization:

.. math::

    SL_S=L_S\Lambda_S

and then the orthogonalization matrix is constructed:

.. math::

    S_{ortho}=L_S\Lambda^{−1/2}_SL^T_S

The SCF iterations requires an initial Fock matrix given as:

.. math::
    
    F_0=H^{core}

The SCF procedure is calulated as the following equations. First the fock matrix is constructed:

.. math::

    F_{n,ij}∣_{n≠0} = H^{core}_{ij}
                    +\sum_{kl}^{AO} D_{n−1,kl} (2(ij||kl)−(ik||jl))

Then the Fock matrix is brought into the orthogonal basis:

.. math::

    F^{\prime}_n=S^T_{ortho} H_{core}S_{ortho}

The F’ is then diagonalized:

.. math::

    F^{\prime}_nC^{\prime}_n=C^{\prime}_n \epsilon_n

The coefficients are transformed back:

.. math::

    C_n=S_{ortho}C^{\prime}_n

A density matrix can be made from the coefficients:

.. math::

    D_{n,ij}=\sum_k^{occ} C_{n,ki}C_{n,kj}

The electronic energy of system can be found as:

.. math::

    E_{n,elec}= \sum_{ij}^{AO} D_{0,ij}(H^{core}_{ij}+F_{n,ij})

The above SCF procedure will is stopped at certain tresholds. The change in energy and the RMSD
of the density matrix can be found as:

.. math::

    \Delta E_n &= E_{n,elec}−E_{n−1,elec}\\
    RMSD_n     &= \sqrt{\sum_{ij} D_{n,ij}−D_{n−1,ij}}

Plain Hartree-Fock solver
    .. literalinclude:: ../data/examples/hf/scf.py
        :caption: /data/examples/hf/scf.py

DIIS Solver
============
Direct Inversion in the Iterative Subspace (DIIS). Makes new :math:`F^{\prime}` guesses bassed on previous guesses.

The error vector in DIIS is given as:

.. math::

    e_i=F_iD_iS−SD_iF_i

It is wanted that the sum of error vectors is zero:

.. math::

    e^{\prime}=\sum_i c_ie_i=0

And now with the requirement that the sum of all c is zero, the following matrix eqution is solwed:

.. math::

    B_{i,j}c_i=b_0

Here:

.. math::

    B_{ij}=tr(e_i⋅e_j)

and,

.. math::

 b_0 =  \begin{pmatrix}
        0\\
        0\\
        \vdots\\
        -1
        \end{pmatrix}


Finally the new F’ is constructed as:

.. math::

    F^{\prime}=\sum_i c_iF_i

DIIS Hartree-Fock solver
    .. literalinclude:: ../data/examples/hf/scf_diis.py
        :caption: /data/examples/hf/scf_diis.py
