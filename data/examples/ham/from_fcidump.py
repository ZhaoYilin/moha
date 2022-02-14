from moha import *
import numpy as np

fcidump = FCIDUMP(pg='d2h').read('./data/N2.STO3G.FCIDUMP')
ham = ChemicalHamiltonian.from_fcidump(fcidump)
Enuc = ham.operators[OperatorNames.Enuc]
Hcore = ham.operators[OperatorNames.Hcore]
print(Enuc)
print(Hcore)
