============
Hartree-Fock
============
The Hartree--Fock method is an uncorrelated mean-field theory that offers a qualitative description of chemical systems. Although Hartree--Fock theory is only qualitatively correct, it forms the basis for more accurate models and becomes the cornerstone of ab initio quantum chemistry.

Hartree--Fock method assumes that the ground state wave function :math:`\vert\Phi_0\rangle` is an antisymmetrized product state with a single--electron orbital basis :math:`\{\phi_i\}`. It is often represented by a single Slater determinant:

.. math::

   \Phi_0(x_1,x_2,\dots,x_N) = \frac{1}{(\sqrt{N!})}
   \begin{vmatrix}
        \phi_1(x_1) & \phi_2(x_1) & \dots & \phi_N(x_1) \\
        \phi_1(x_2) & \phi_2(x_2) & \dots & \phi_N(x_2) \\
        \vdots      & \vdots      & \dots & \vdots \\
        \phi_1(x_N) & \phi_2(x_N) & \dots & \phi_N(x_N) \\
   \end{vmatrix}

To optimize the Hartree--Fock wave function, we invoke the variation principle to minimize the expectation value of the energy:

.. math::
    
   E_{HF} = \mathop{min}\limits_{\{\phi_i\}}^{} \langle\Phi_0\vert\hat{H}\vert\Phi_0\rangle

The minimization of the expected energy with a given basis set leads to the Roothan--Hartree--Fock equation:   

.. math::

    \mathbf{F(C)C} = \mathbf{SC\epsilon}

where :math:`\mathbf{S}` is the overlap matrix with atomic basis set, :math:`\mathbf{C}` is the matrix of molecular orbital coefficients, :math:`\mathbf{\epsilon}` is a diagonal matrix whose elements represent the eigenenergies, and the :math:`\mathbf{F}` is the Fock matrix, which depends on the coefficients matrix :math:`\mathbf{C}`. 

The Fock matrix :math:`\mathbf{F}` is:

.. math::

   \mathbf{F(C)} = \mathbf{H} + \mathbf{G(C)}

:math:`\mathbf{H}` is the core--Hamiltonian matrix, and it is the sum of kinetic energy integral :math:`\mathbf{T}` and electron--nuclear attraction integral :math:`\mathbf{V}`. Both are one--electron integral and invariance during the self--consistent field optimization.   

.. math::

   \mathbf{H} = \mathbf{T} + \mathbf{V}

:math:`\mathbf{G}` is Hartree--Fock potential matrix, and it is the sum of coulomb matrix :math:`\mathbf{J}` and the exchange matrix :math:`\mathbf{K}`. Both of them are two--electron integral matrix and is a function of :math:`\mathbf{C}`.

.. math::

   \mathbf{G(C)} = \mathbf{J(C)} + \mathbf{K(C)}

Hartree--Fock energy is an upper bound to the exact energy using the variational theorem. It typically gives errors of :math:`0.5\%` in the energy, :math:`1\%` in molecular geometry, and :math:`5−10\%` in other properties. The correlation energy is the difference between the Hartree--Fock limit energy and a system's exact non--relativistic energy.

.. math::

    E_{HF} &\geq E_{exact}\\
    E_{corr} &= E_{exact} - E_{HF}

Plain Solver
============
Since the Hartree--Fock--Roothaan equations are nonlinear, we will solve them with an iterative procedure. Such a procedure is often referred to as a self--consistent field (SCF) algorithm, and here we introduce the most naive one.

* **Step 1: Build The Integral**

The pre--computed quantities for the SCF algorithm include:

    * The nuclear repulsion energy :math:`E_{nuc}`.

    * One--electron integrals overlap :math:`\mathbf{S}`, kinetic energy :math:`\mathbf{T}`, and nuclear attraction :math:`\mathbf{V}`.

    * Two--electron integral, the electron repulsion :math:`(\mu\nu\vert\lambda\sigma)`.

Please check out the :ref:`Hamiltonian <Hamiltonian>` section for further information to learn more about the pre--computed quantities.    

* **Step 2: Build the Orthogonalization Matrix**

Hartree--Fock algorithm must be optimized with the constraint that the orbitals are orthonormal, ensuring that the Hartree--Fock wave function is normalized.  

.. math::

   \langle \phi_i \vert \phi_j \rangle = \delta_{ij}

The atomic orbitals are, in general, not orthonormal. Hence we bring the orthogonalization matrix into the algorithm to make,

.. math::

   \mathbf{X^{\dagger}SX} = \mathbf{I}

where :math:`\mathbf{X}` is the orthogonalization matrix, and :math:`\mathbf{I}` is the density matrix.

The orthogonalization matrix :math:`\mathbf{X}` can be constructed from the overlap matrix :math:`\mathbf{S}` by diagonalization:

.. math::

   \mathbf{s} = \mathbf{USU}^{\dagger}

where :math:`\mathbf{s}` are eigenvalues of the overlap matrix :math:`\mathbf{S}`, further :math:`\mathbf{U}` is a matrix of the eigenvector of the overlap matrix :math:`\mathbf{S}`.

The construction of orthogonalization matrix :math:`\mathbf{X}` is not unitary, and there are two ways of orthogonalization are most frequently used.

The first one, called symmetric orthogonalization:

.. math::

    \mathbf{X} = \mathbf{S}^{-1/2} = \mathbf{Us}^{-1/2}\mathbf{U}^{\dagger}

The second way is canonical orthogonalization:

.. math::

   \mathbf{X} = \mathbf{Us}^{-1/2}

* **Step 3: Set the Initial Guess**

The SCF iterations require an initial guess of the coefficients matrix. As the naive solver, we give zero matrices as the initial guess by default,  

.. math::

   \mathbf{C}_0 = \mathbf{0}

The naive initial guess of coeffiences matrix :math:`\mathbf{C}_0` results in zero matrix as :math:`\mathbf{D}_0` and core--Hamiltonian as initial Fock matrix:

.. math::

   \mathbf{D}_0 = \mathbf{0}

.. math::
    
    \mathbf{F}_0 = \mathbf{H}_{core} + \mathbf{G}(\mathbf{0}) = \mathbf{H}_{core}

* **Step 4: Calculate density matrix from Fock matrix**

To calculate the new density matrix :math:`\mathbf{D}_i` (i denotes the current iteration number), first bring the Fock matrix :math:`\mathbf{F}_{i-1}` into the orthogonal basis, and we get the transformed Fock matrix :math:`\mathbf{F}^{\prime}_{i-1}`:

.. math::

    \mathbf{F}^{\prime}_{i-1} = \mathbf{X}^T \mathbf{F}_{i-1} \mathbf{X}

The transformed Fock matrix :math:`\mathbf{F}^{\prime}_{i-1}` is then diagonalized:

.. math::

    \mathbf{F}^{\prime}_{i-1} \mathbf{C}^{\prime}_{i} = \mathbf{C}^{\prime}_{i} \mathbf{\epsilon}_i

Back-transform:    

.. math::

    \mathbf{C}_i = \mathbf{X}\mathbf{C}^{\prime}_i

Calculate the density matrix:

.. math::

    \mathbf{D}_i = \mathbf{C}_{i} \mathbf{n} \mathbf{C}^{\dagger}_{i} 

where :math:`\mathbf{n}` is the electron occupation matrix.

* **Step 5: Calculate Fock matrix from density matrix**

With the new density matrix :math:`\mathbf{D}_i`, we can compute the current iteration Fock matrix :math:`\mathbf{F}_i` by:

.. math::

    \mathbf{F}_{i,\mu\nu} = \mathbf{h}_{\mu\nu} 
    + \sum_{\sigma,\tau} (2\langle\mu\sigma\vert\vert\nu\tau\rangle−\langle\mu\sigma\vert\vert\nu\tau\rangle)\mathbf{D}_{i,\sigma\tau}

* **Step 6: Self-consistency check converagence**

Repeat steps 4 and 5, compares the change of electronic energy and-or the root-mean-squared difference (RMSD) of the density matrix between two iterations. Furthermore, the above SCF procedure will stop at certain thresholds.

The change in the electronic energy and RMSD of the density matrix can be found as:

.. math::

   \Delta E_{i,elec} &= E_{i,elec} - E_{i-1,elec} \leqslant \sigma_1\\
   RMSD_i &= \sqrt{\sum_{\mu\nu} D_{i,\mu\nu}−D_{i−1,\mu\nu}} \leqslant \sigma_2

* **Step 7: Calculate energy**

The electronic energy of system can be found as:

.. math::

    E_{elec}= \sum_{\mu\nu}^{AO} D_{\mu\nu}(h_{\mu\nu}+F_{\mu\nu})

The total energy is the sum of electronic energy and the nuclear repulsion energy:    
    
.. math::

    E_{tot}= E_{elec} + E_{nuc} 

Plain Hartree-Fock solver
    .. literalinclude:: ../../data/examples/hf/scf.py
        :caption: /data/examples/hf/scf.py

DIIS Solver
============
The plain Hartree--Fock solver updates orbitals from the Fock operator in last iteration. It may fail to converge or converge slowly. The Direct Inversion in the Iterative Subspace (DIIS) solver shares the same algorithm with the plain solver but generates orbitals from an averaged Fock operator of preceding iterations hence converging fast.
   
* **Stpe 1: Claculate the Error Matrix in Each Iteration**

The error vector in DIIS is given as:

.. math::

    e_i=F_iD_iS−SD_iF_i

It is wanted that the sum of error vectors is zero:

.. math::

    e^{\prime}=\sum_i c_ie_i=0

where :math:`c_i` is coefficient of the normalize constraint

.. math::

   \sum_i c_i = 1

* **Stpe 2: Build the B matrix and Solve the Linear Equations**

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


* **Stpe 3: Compute the New Fock Matrix**

With a set of Fock matrices :math:`\{\mathbf{F}_i\}`, the new Fock matrix :math:`\mathbf{F}^{\prime}` is constructed as:

.. math::

    F^{\prime}=\sum_i c_iF_i



DIIS Hartree-Fock solver
    .. literalinclude:: ../../data/examples/hf/scf_diis.py
        :caption: /data/examples/hf/scf_diis.py
