Molecular System
################
To begin a calculation with MoHa, the first step is to build a Hamiltonian of a 
system, either molecular system or model system. In most cases, we work with a molecular
system.

In terms of second quantisation operators, a general Hamiltonian can be written
as 

.. math::
        H = - \sum_{ij} t_{ij}\hat{c}^{\dagger}_{i}\hat{c}_{j} + \frac{1}{2} \sum_{ijkl}
        V_{ijkl}\hat{c}^{\dagger}_{i}\hat{c}^{\dagger}_{k}\hat{c}_{l}\hat{c}_{j}

The construction of molecular Hamiltonian usually set up in three steps. 

- First, construct a molecular geometry. 
- Second, generate a basis set for the molecular.
- Finally, compute all kinds of one body terms and two body terms with that basis
  to define a Hamiltonian.

 
Molecular Geometry
==================
Molecule is a system consist with nucleus and electrons. For quantum chemistry
calculation in MoHa, we will always used the Born-Oppenheimer approximation, which assume
that the motion of atomic nuclei and electrons in a molecule can be separated

.. math::
    \Psi_{molecule} = \psi_{electronic} \otimes \psi_{nuclear}

The module ``molecule`` in MoHa actually only contains imformation of the
nuclear. To build a water molecule with MoHa,  we can specify the molecular object by
load the molecular geometry from .xyz file.

.. code-block:: text

    3
    water
    8   0.0000000000	-0.143225816552		0.0000000
    1   1.638036840407	1.136548822547		-0.000000000000
    1   -1.638036840407	1.136548822547		-0.000000000000

.. code-block:: python

    mol = IOSystem.from_file('h2o.xyz')

.. figure:: ./pictures/water.png
    :scale: 40 %
    :align: center

Basis Set
=========
MoHa supports basis sets consisting of generally contracted Cartesian Gaussian
functions. MoHa is using the same basis set format as NWChem, and the basis sets can be downloaded from the EMSL webpage (https://bse.pnl.gov/bse/portal).

STO-3G EMSL basis set of hydrogen and oxygen:

.. code-block:: text

    #  STO-3G  EMSL  Basis Set Exchange Library  8/21/18 3:24 AM
    # Elements                             References
    # --------                             ----------
    #  H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 2657 (1969).
    # Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart, J.A. Pople,
    #          J. Chem. Phys.  2769 (1970).
    # K,Ca - : W.J. Pietro, B.A. Levy, W.J. Hehre and R.F. Stewart,
    # Ga - Kr: J. Am. Chem. Soc. 19, 2225 (1980).
    # Sc - Zn: W.J. Pietro and W.J. Hehre, J. Comp. Chem. 4, 241 (1983) + Gaussian.
    #  Y - Cd: W.J. Pietro and W.J. Hehre, J. Comp. Chem. 4, 241 (1983). + Gaussian
    #   


    BASIS "ao basis" PRINT
    #BASIS SET: (3s) -> [1s]
    H    S
        3.42525091             0.15432897       
        0.62391373             0.53532814       
        0.16885540             0.44463454       
    #BASIS SET: (6s,3p) -> [2s,1p]
    O    S
        130.7093200              0.15432897       
        23.8088610              0.53532814       
        6.4436083              0.44463454       
    O    SP
        5.0331513             -0.09996723             0.15591627       
        1.1695961              0.39951283             0.60768372       
        0.3803890              0.70011547             0.39195739       
    END

.. code-block:: python

    mol,orb = IOSystem.from_file('h2o.xyz','sto-3g.nwchem')

Hamiltonian
==========
The core mechanical quantities of a chemistry system is the Hamiltonian. Hamiltonian operator
should include the kinetic energy and potential energy terms of all atomic nuclei and all
electrons. It is generally assumed that the molecule is in a vacuum and adiabatic state
in isolation. At this time, the interaction potential energy between the nucleus and the
electron in the molecule is only related to distance from each other and time independent. 
Its expression is:

.. math::
    \hat{H}= &-\sum^N_{i=1}\frac{\hbar^2}{2m_i}{\nabla}_i^2
        -\sum^N_{i=1}\frac{\hbar^2}{2M_\alpha}{\nabla}_\alpha^2
        - \sum^N_{i=1}\sum^M_{\alpha=1} \frac{Z_\alpha e^2}{\textbf{r}_{i\alpha}}\\
        &+ \sum^N_{\alpha=1}\sum^M_{\beta=1} \frac{Z_\alpha Z_\beta e^2}{\textbf{R}_{\alpha\beta}}
        +\sum^N_{i=1}\sum^N_{j>i} \frac{e^2}{\textbf{r}_{ij}}

The formula contains five terms:

Kinetic energy of electrons.

.. math::
  		\hat{T}_e = -\sum^N_{i=1}\frac{\hbar^2}{2m_i}\boldsymbol{\nabla}_i^2

Kinetic energy of atomic nuclei.

.. math::
  		\hat{T}_n = -\sum^N_{i=1}\frac{\hbar^2}{2M_\alpha}{\nabla}_\alpha^2 

Nuclear attraction.

.. math::
  		\hat{V}_{en} = -\sum^N_{i=1}\sum^M_{\alpha=1} \frac{Z_\alpha e^2}{\textbf{r}_{i\alpha}}  

Repulsive between nuclei.

.. math::
  		\hat{V}_{nn} = \sum^N_{\alpha=1}\sum^M_{\beta=1} \frac{Z_\alpha Z_\beta e^2}{\textbf{R}_{\alpha\beta}} 

Repulsive between electrons.

.. math::
  		\hat{V}_{ee} = \sum^N_{i=1}\sum^N_{j>i} \frac{e^2}{\textbf{r}_{ij}}

:math:`m_i` is the mass of electron. :math:`M_\alpha` and :math:`Z_\alpha` refer to the mass and charge of atomic nucleus. 
:math:`R_{\alpha\beta}`, :math:`r_{i\alpha}` and :math:`r_{ij}` is the distance between two nucleus, atomic nuclei 
and electron and two electrons respectively. The explicit representation of Laplacian operator is:

.. math::
	\boldsymbol{\nabla}^2 = \frac{\partial^2}{\partial x^2} +\frac{\partial^2}{\partial y^2} 
	+ \frac{\partial^2}{\partial z^2}


To build a Hamiltonian object, MoHa can load the molecular geometry and and basis from file
format.

.. code-block:: python

    mol,orbs = IOSystem.from_file('h2o.xyz','sto-3g.nwchem')
    ham = Hamiltonian.build(mol,orbs)

Hamiltonian object has attributes of different operators use the following
conventions for variable names. The following are defined by setting up the
Hamiltonian by default:

* ``ham.operators['nuclear_repulsion']``
    :class:`.ZeroElectronOperator` object with the nuclear repulsion energy integrals.
* ``ham.operators['overlap']``
    :class:`.OneElectronOperator` object with the overlap integrals.
* ``ham.operators['kinetic']``
    :class:`.OneElectronOperator` object with the kinetic energy integrals.
* ``ham.operators['nuclear_attraction']``
    :class:`.OneElectronOperator` object with the nuclear attraction integrals.
* ``ham.operators['electron_repulsion']``
    :class:`.TwoElectronOperator` object with the electron repulsion integrals.

They offer the key ingredient for the following calculations.
