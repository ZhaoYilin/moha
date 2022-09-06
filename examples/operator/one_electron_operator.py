from moha import *
import numpy as np

mol,orbs = IOSystem.from_file('../data/water.xyz','sto-3g.nwchem')
S = OverlapOperator.build(orbs)
T = KineticOperator.build(orbs)
V = NuclearAttractionOperator.build(mol,orbs)
print(S)
print(T)
print(V)
