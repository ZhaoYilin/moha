import numpy as np
from moha.io.log import log,timer
from moha.system.operator import MultipoleMomentOperator

class MultipoleMoment(object):
    """Multipole momnet solver.

    Attributes
    ----------
    mol
        Chemical molecule.

    basis_set
        Basis set of the molecule

    wfn
        Hartree Fock wavefunction.

    Methods
    -------
    __init__(self,mol,basis_set,ham,wfn)
        Initialize the solver.

    assign_molecule(self,mol)
        Assign chemical molecule to the solver.

    assign_basis_set(self,basis_set)
        Assign basis set to the solver.

    assign_wavefunction(self,wfn)
        Assign Hartree Fock wavefunction to the solver.
    """
    def __init__(self,mol,basis_set,wfn):
        """Initialize the solver.

        Attributes
        ----------
        mol
            Chemical molecule.

        basis_set
            Basis set of the molecule

        wfn
            Hartree Fock wavefunction.
        """
        self.assign_molecule(mol)
        self.assign_basis_set(basis_set)
        self.assign_wavefunction(wfn)

    def assign_molecule(self,mol):
        """Assign chemical molecule to the solver.

        Attributes
        ----------
        mol
            Chemical molecule.
        """
        self.mol = mol

    def assign_basis_set(self,basis_set):
        """Assign basis set to the solver.

        Attributes
        ----------
        basis_set
            Basis set of the molecule.
        """
        self.basis_set = basis_set

    def assign_wavefunction(self,wfn):
        """Assign Hartree Fock wavefunction to the solver.

        Attributes
        ----------
        wfn
            Hartree Fock wavefunction.
        """
        self.wfn = wfn

    @timer.with_section('Multipole')
    def kernel(self,lmns,coordinate = np.array([0.,0.,0.])):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            Multipole moment results.
        """
        log.hline()
        log('Multipole moment Section'.format())
        log.hline()
        
        D = self.wfn.density_matrix
        mms = []

        for lmn in lmns:
            mm_matrix = MultipoleMomentOperator.build(self.basis_set,lmn,coordinate).integral
            mm = np.dot(D,mm_matrix)
            mm = np.trace(mm,axis1=0,axis2=1)
            mms.append(mm)

        log.hline()
        log('Multipole Moment'.format())
        log('{}'.format(mms))
        log.hline()

        results = {
        "success": True,
        "mms": mms
        }
        return results

