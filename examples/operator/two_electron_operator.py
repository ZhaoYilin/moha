from moha import *
import numpy as np

mol,orbs = IOSystem.from_file('../data/water.xyz','sto-3g.nwchem')
Eri = ElectronRepulsionOperator.build(orbs)
print(Eri)
