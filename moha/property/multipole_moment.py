import numpy as np
from moha.system.operator import OneElectronOperator
class OneElectronProperty(object):
    def __init__(self):
        pass
    @classmethod
    def multipole_moment(cls,mol,orbs,wavefunction,lmns,coordinate = np.array([0.,0.,0.])):
        D = wavefunction.density
        mms = []
        for lmn in lmns:
            mm_matrix = OneElectronOperator.build('multipole_moment',mol,orbs,lmn,coordinate).value
            mm = np.dot(D,mm_matrix)
            mm = np.trace(mm,axis1=0,axis2=1)
            mms.append(mm)
        return mms
