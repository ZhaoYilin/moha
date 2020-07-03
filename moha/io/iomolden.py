import sys
import numpy as np
import moha

class Section(object):
    """
    store the label of section
    """
    def __init__(self):
        pass



class Molden(object):
    """
    """
    def __init__(self,filename,mol,basis,wf):
        """
        """
        self.filename = filename
        self.mol = mol
        self.basis = basis
        self.wf = wf

    def dump_molden(self):
        """write molden input file
        
        **Arguments:**

        filename
            File name of molden input file

        data
        """
        with open(self.filename,'w') as f:
            f.write("[Molden Format]")
            f.write("[Title]")
            f.write("[Molden Format]")
            f.write()
            #self.write_atoms()
            #self.write_GTO()
            #self.wrtie_MO()

    def write_atoms(self):
        with open(filename,'w') as f:
            f.write("[Atoms] AU")
            for i in range():
    def write_GTO(self):
        pass
    def write_MO(self):
        pass



