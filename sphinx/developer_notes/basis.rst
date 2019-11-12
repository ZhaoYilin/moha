Basis Set
=========
MoHa supports basis sets consisting of generally contracted Cartesian Gaussian
functions. MoHa is using the same basis set format as NWChem, and the basis sets can be downloaded from the EMSL webpage (https://bse.pnl.gov/bse/portal).

The basis object is a list of list with the follwing structure

.. math::

    [Idx x y z #CONTR CONTR Atomidx]

Idx is the AO index. Atomidx, is the index of the associated atom. #CONTR contains the number of primitive functions in the contracted.

CONTR contains all the information about the primitive functions and have the form:

.. math::

    [N ζ c l m n]

Inside the integral call, the basisset file is reconstructed into three different arrays, containing the basisset information. The first one is basisidx that have the following form:

.. math::

    [#primitives\ \ \ loopstartidx]

It thus contains the number of primitives in each basisfunction, and what start index it have for loop inside the integral code.

The second array is basisint, that have the following forms:

.. math::

    [l m n]
    [l m n atomidx]

The first one is for regular integrals and the second one is for derivatives. Both contains all the angular momentum quantum numbers, and the derivative also contains the atom index (used in derivative of VNe).

The last array is basisfloat and have the following forms:

.. math::
    [N ζ c x y z]
    [N ζ c x y z N_{x,+} N_{x,−} N_{y,+} N_{y,−} N_{z,+} N_{z,−}]

basisfloat contains the normalization constants, Gaussian exponent and prefacor and the coordinates of the atoms. The second one is again for the derivatives, it contains normalization constants of the differentiated primitives.

