======
System
======
The ``system`` module is the core of MoHa, and almost every MoHa calculation starts from this module. To begin a MoHa calculation, first build the Hamiltonian of a system. In terms of second quantization operators, a time-independent non-relativistic Hamiltonian gives:

.. math::
        H = - \sum_{ij} t_{ij}\hat{c}^{\dagger}_{i}\hat{c}_{j} + \frac{1}{2} \sum_{ijkl}
        V_{ijkl}\hat{c}^{\dagger}_{i}\hat{c}^{\dagger}_{k}\hat{c}_{l}\hat{c}_{j}

The construction of the molecular Hamiltonian is set in three steps.

- Construct a molecular geometry.
- Generate a basis set for the molecular.
- Generate all kinds of operator terms with molecule and basis set to define a Hamiltonian.
 
Molecular Geometry
==================
A molecule is a system consisting of a nucleus and electrons. For quantum chemistry calculation in MoHa, we will always use the Born-Oppenheimer approximation, which assumes that the motion of atomic nuclei and electrons in a molecule can be separated.

.. math::

    \Psi_{molecule} = \psi_{electronic} \otimes \psi_{nuclear}

The class :class:`Molecule` in MoHa only contains information about the nuclear. Therefore, to build a water molecule with MoHa, we can specify the molecular object by loading the molecular geometry from xyz file.    

The XYZ file format is a chemical file format. There is no formal standard and several variations exist, but a typical XYZ format specifies the molecule geometry by giving the number of atoms with Cartesian coordinates that will be read on the first line, a comment on the second, and the lines of atomic coordinates in the following lines.

* Format

The formatting of xyz is:

.. code-block:: text

    <number of atoms>
    comment line
    <element>   <X>    <Y>    <Z>
    \dots

* Example

The water molecule in the xyz format:

.. code-block:: text

    3
    water
    8   0.0000000000	-0.143225816552		0.0000000
    1   1.638036840407	1.136548822547		-0.000000000000
    1   -1.638036840407	1.136548822547		-0.000000000000

.. code-block:: python

    import moha

    mol = moha.io.iogeometry.load_xyz('h2o.xyz')

The first line in the xyz file is the number of atoms in the molecule, while the second line gives the name of the molecule, remaining lines give the atomic number and coordinate of each atom. In quantum chemistry, the atomic unit system is generally used, here Bohr radius (:math:`a_0` = 0.0529177nm) is taken as the unit by default.

Basis Set
=========
MoHa supports Cartesian Gaussian type orbtial(GTO), it is an atomic orbital used in linear combinations forming molecular orbitals. Cartesian GTOs are defiend by an angular par which is homogeneous polynomial in the compoents x, y, and z of the position vector :math:`\mathbf{r}`. That is,

.. math::
    
    G_{ijk}(r_A,\alpha) = x^i_A y^j_A z^k_A e^{-\alpha r^2_A}

here :math:`r_A = r - A` is the atom-centered function :math:`\alpha > 0` is the real orbital exponent and :math:`i+j+k = n`.

MoHa use the same basis set input file as NWChem, all the available basis set file can be found in the :file:`moha/data/examples` directory. If you need basis set more than MoHa offered, you can downloaded it from the EMSL webpage (https://bse.pnl.gov/bse/portal).

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

    from moha.io.iosystem import IOSystem

    mol,orb = IOSystem.from_file('h2o.xyz','sto-3g.nwchem')

.. _Hamiltonian:
    
Hamiltonian
===========
The core mechanical quantities of a chemistry system is the Hamiltonian. Hamiltonian operator should include the kinetic energy and potential energy terms of all atomic nuclei and all electrons. It is generally assumed that the molecule is in a vacuum and adiabatic state in isolation. At this time, the interaction potential energy between the nucleus and the electron in the molecule is only related to distance from each other and time independent. Its expression is:

.. math::
    \hat{H}= &-\sum^N_{i=1}\frac{\hbar^2}{2m_i}{\nabla}_i^2
    - \sum^N_{i=1}\sum^M_{\alpha=1} \frac{Z_\alpha e^2}{\textbf{r}_{i\alpha}}\\
    &+\sum^N_{i=1}\sum^N_{j>i} \frac{e^2}{\textbf{r}_{ij}}
    +\sum^N_{\alpha=1}\sum^M_{\beta=1} \frac{Z_\alpha Z_\beta e^2}{\textbf{R}_{\alpha\beta}}

The formula contains four terms:

* Kinetic energy of electrons.

.. math::

    \hat{T}_e = -\sum^N_{i=1}\frac{\hbar^2}{2m_i}\boldsymbol{\nabla}_i^2

* Nuclear attraction.

.. math::
  		
    \hat{V}_{en} = -\sum^N_{i=1}\sum^M_{\alpha=1} \frac{Z_\alpha e^2}{\textbf{r}_{i\alpha}}  

* Repulsive between electrons.

.. math::
  		
    \hat{V}_{ee} = \sum^N_{i=1}\sum^N_{j>i} \frac{e^2}{\textbf{r}_{ij}}

* Repulsive between nuclei.

.. math::
  		
    \hat{V}_{nn} = \sum^N_{\alpha=1}\sum^M_{\beta=1} \frac{Z_\alpha Z_\beta e^2}{\textbf{R}_{\alpha\beta}} 


:math:`m_i` is the mass of electron. :math:`M_\alpha` and :math:`Z_\alpha` refer to the mass and charge of atomic nucleus. :math:`R_{\alpha\beta}`, :math:`r_{i\alpha}` and :math:`r_{ij}` is the distance between two nucleus, atomic nuclei and electron and two electrons respectively. The explicit representation of Laplacian operator is:

.. math::
	\boldsymbol{\nabla}^2 = \frac{\partial^2}{\partial x^2} +\frac{\partial^2}{\partial y^2} 
	+ \frac{\partial^2}{\partial z^2}


To build a Hamiltonian object, MoHa need both molecular geometry and basis object.

.. code-block:: python

    from moha.io.iosystem import IOSystem
    from moha.system.hamiltonian.chemical_hamiltonian import ChemicalHamiltonian

    mol,orbs = IOSystem.from_file('h2o.xyz','sto-3g.nwchem')
    ham = ChemicalHamiltonian.build(mol,orbs)

Hamiltonian object has attributes of different operators use the following conventions for variable names. The following are defined by setting up the Hamiltonian by default:

.. code-block:: python
    
    
    from moha.io.iosystem import IOSystem
    from moha.system.hamiltonian.chemical_hamiltonian import ChemicalHamiltonian
    from moha.system.operator.base import OperatorNames

    mol,orbs = IOSystem.from_file('h2o.xyz','sto-3g.nwchem')
    ham = ChemicalHamiltonian.build(mol,orbs)

    nuclear_energy = ham.operators[OperatorNames.Enuc]
    overlap = ham.operators[OperatorNames.S]
    kinetic = ham.operators[OperatorNames.T]
    nuclear_attraction = ham.operators[OperatorNames.V]
    electron_repuslion = ham.operators[OperatorNames.Eri]

They offer the key ingredient for the following calculations.
