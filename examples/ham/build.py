from moha import *

import numpy as np

geo = [[8,   0.000000000000,  -0.143225816552,   0.000000000000],
    [1,   1.638036840407,   1.136548822547,  -0.000000000000],
    [1,  -1.638036840407,   1.136548822547,  -0.000000000000]]

mol = Molecule.build(geo)
orb = BasisSet.build(mol,'sto-3g.nwchem')

ham = ChemicalHamiltonian.build(mol,orb)
