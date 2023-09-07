import numpy as np
from moha.io.log import log,timer
from moha.hamiltonian.operator.multipole_moment import MultipoleMoment

class Moment(object):
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

    @timer.with_section('Moment')
    def multipole(self,lmns,coordinate = np.array([0.,0.,0.])):
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
            mm_matrix = MultipoleMoment.build(self.basis_set,lmn,coordinate)
            mm = -np.dot(D,mm_matrix)
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

    @timer.with_section('Moment')
    def dipole(self,coordinate = np.array([0.,0.,0.])):
        """Kernel of the solver.

        Returns
        -------
        results : dict
            Multipole moment results.
        """
        log.hline(char='=')
        log.blank()
        log('Dipole moment Section'.format())

        lmns = [[1,0,0],[0,1,0],[0,0,1]]
        D = self.wfn.density_matrix
        mms = []

        for lmn in lmns:
            mm_matrix = MultipoleMoment.build(self.basis_set,lmn,coordinate)
            mm = -np.dot(D,mm_matrix)
            mm = np.trace(mm,axis1=0,axis2=1)
            nuclear_mm = 0
            for atom in self.mol:
                x,y,z = atom.coordinate
                l,m,n = lmn
                nuclear_mm += atom.number*x**l*y**m*z**n
            mm = mm + nuclear_mm
            mms.append(mm)

        log.hline()
        log('Multipole Moment'.format())
        log('{0:<20s}{1:<20.6f}'.format('X',mms[0]))
        log('{0:<20s}{1:<20.6f}'.format('Y',mms[1]))
        log('{0:<20s}{1:<20.6f}'.format('Z',mms[2]))
        log('{0:<20s}{1:<20.6f}'.format('Total',sum(mms)))
        log.blank()

        results = {
        "success": True,
        "mms": mms
        }

        return results

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
