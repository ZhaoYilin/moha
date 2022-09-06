import numpy as np
import re

class FCIDUMP(object):
    """An interface between python and FCIDUMP format file.

    Attributes
    ----------
    n_orb: int
        The number of orbitals.

    n_elec: int
        The number of electrons.

    ms2: int
        The spin polarization of the system.

    orb_sym: list
        A list containing the symmetry label of each orbital.

    isym: int
        Total symmetry of the wavefunction.

    e_nuc: float
        Scalar representation of the nuclear repulsion energy.

    epsilon: list
        Vector representation of orbital energies.

    h1e: np.ndarray
        Matrix representation of one electron operator.

    g2e: np.ndarray
        Tensor representation of two electron operator.

    Methods
    -------
    __init__(self, n_orb=0, n_elec=0, ms2=0, orb_sym=None,
                 isym=1, e_nuc=0, epsilon=None, h1e=None, g2e=None)
        Initialize the instance.

    load(self, filename)
        Load fci information from FCIDUMP file.

    dump(self, filename, tol=1E-13)
        Dump fci information to FCIDUMP file.
    """
    def __init__(self, n_orb=0, n_elec=0, ms2=0, orb_sym=None,
                 isym=1, e_nuc=0, epsilon=None, h1e=None, g2e=None):
        """Initialize the instance.

        Parameters
        ----------
        n_orb: int
            The number of orbitals.

        n_elec: int
            The number of electrons.

        ms2: int
            The spin polarization of the system.

        orb_sym: list
            A list containing the symmetry label of each orbital.

        isym: int
            Total symmetry of the wavefunction.

        e_nuc: float
            Scalar representation of the nuclear repulsion energy.

        epsilon: list
            Vector representation of orbital energies.

        h1e: np.ndarray
            Matrix representation of one electron operator.

        g2e: np.ndarray
            Tensor representation of two electron operator.
        """
        self.n_orb = n_orb
        self.n_elec = n_elec
        self.ms2 = ms2

        self.orb_sym = [0] * self.n_orb if orb_sym is None else orb_sym
        self.isym = isym

        self.epsilon = epsilon
        self.e_nuc = e_nuc
        self.h1e = h1e
        self.g2e = g2e

    def load(self, filename):
        """Load fci information from FCIDUMP file.

        Parameters
        ----------
        filename : str
            Name of the fcidump file.
    
        Return
        ------
        fcidump : FCIDUMP instance
            A FCIDUMP instance.
        """
        with open(filename, 'r') as handle:
            content = handle.read().lower()
            if '/' in content:
                pars, ints = content.split('/')
            elif '&end' in content:
                pars, ints = content.split('&end')
            handle.close()

        #Assign parameters
        par_dict = {}
        pars = pars.replace('&fci','').replace(' ','').replace('\n','').replace(',,',',')
        for par in re.split(',(?=[a-zA-Z])', pars):
            key, val = par.split('=')
            if key in ('norb', 'nelec', 'ms2', 'isym'):
                par_dict[key] = int(val.replace(',', ''))
            elif key in ('orbsym'):
                par_dict[key] = [int(x) for x in val.replace(',', ' ').split()]
            else:
                par_dict[key] = val
        self.n_orb = int(par_dict.get('norb'))
        self.n_elec = int(par_dict.get('nelec', 0))
        self.twos = int(par_dict.get('ms2', 0))
        self.isym = int(par_dict.get('isym', 1))
        self.orb_sym = par_dict.get('orbsym')

        #Assign integrals
        data = []
        dtype = float
        for line in ints.split('\n'):
            line = line.strip()
            if len(line) == 0 or line.strip()[0] == '!':
                continue
            line = line.split()
            assert len(line) == 5 or len(line) == 6
            if len(line) == 5:
                value = float(line[0])
                i, j, k, l = [int(x) for x in line[1:]]
            else:
                dtype = np.complex128
                value = float(line[0]) + float(line[1]) * 1j
                i, j, k, l = [int(x) for x in line[2:]]
            data.append((value,i, j, k, l))

        self.epsilon = np.zeros((self.n_orb), dtype=dtype)
        self.h1e = np.zeros((self.n_orb, self.n_orb), dtype=dtype)
        self.g2e = np.zeros((self.n_orb, self.n_orb,
                             self.n_orb, self.n_orb), dtype=dtype)
        for value,i, j, k, l in data:
            if i+j+k+l == 0:
                self.e_nuc = value
            elif j+k+l == 0:
                self.epsilon[i-1] = value
            elif k+l == 0:
                self.h1e[i - 1, j - 1] = self.h1e[j - 1, i - 1] = value
            else:
                self.g2e[i - 1, j - 1, k - 1, l - 1] = value
                self.g2e[j - 1, i - 1, k - 1, l - 1] = value
                self.g2e[j - 1, i - 1, l - 1, k - 1] = value
                self.g2e[i - 1, j - 1, l - 1, k - 1] = value
                self.g2e[k - 1, l - 1, i - 1, j - 1] = value
                self.g2e[k - 1, l - 1, j - 1, i - 1] = value
                self.g2e[l - 1, k - 1, j - 1, i - 1] = value
                self.g2e[l - 1, k - 1, i - 1, j - 1] = value
        return self

    def dump(self, filename, tol=1E-13):
        """Dump fci information to FCIDUMP file.

        Parameters
        ----------
        filename : str
            Name of the fcidump file.

        tol : float 
            Threshold for terms written into file.
        """
        with open(filename, 'w') as handle:
            #Assign parameters
            handle.write(' &FCI NORB=%4d,NELEC=%2d,MS2=%d,\n' %
                       (self.n_orb, self.n_elec, self.ms2))
            if self.orb_sym is not None and len(self.orb_sym) > 0:
                handle.write('  ORBSYM=%s,\n' % ','.join(map(str,self.orb_sym)))
            else:
                handle.write('  ORBSYM=%s\n' % ('1,' * self.n_orb))
            handle.write('  ISYM=%d,\n' % self.isym)
            handle.write(' &END\n')

            #Assign integrals
            output_format = '%20.16f%4d%4d%4d%4d\n'
            npair = self.n_orb * (self.n_orb + 1) // 2

            def write_g2e(handle, g2e):
                assert g2e.ndim in [1, 4]
                if g2e.ndim == 4:
                    # general
                    assert(g2e.size == self.n_orb ** 4)
                    for i in range(self.n_orb):
                        for j in range(self.n_orb):
                            for k in range(self.n_orb):
                                for l in range(self.n_orb):
                                    if abs(g2e[i, j, k, l]) > tol:
                                        handle.write(output_format % (
                                            g2e[i, j, k, l], i + 1, j + 1, k + 1, l + 1))
                else:
                    # 8-fold symmetry
                    assert g2e.size == npair * (npair + 1) // 2
                    ij = 0
                    ijkl = 0
                    for i in range(self.n_orb):
                        for j in range(0, i + 1):
                            kl = 0
                            for k in range(0, i + 1):
                                for l in range(0, k + 1):
                                    if ij >= kl:
                                        if abs(g2e[ijkl]) > tol:
                                            handle.write(output_format % (
                                                g2e[ijkl], i + 1, j + 1, k + 1, l + 1))
                                        ijkl += 1
                                    kl += 1
                            ij += 1

            def write_h1e(handle, h1e):
                for i in range(self.n_orb):
                    for j in range(0, i+1):
                    #for j in range(0, self.n_orb):
                        if abs(h1e[i, j]) > tol:
                            handle.write(output_format %
                                (h1e[i, j], i + 1, j + 1, 0, 0))

            def write_epsilon(handle, epsilon):
                for i in range(self.n_orb):
                    if abs(epsilon[i]) > tol:
                        handle.write(output_format %
                            (epsilon[i], i + 1, 0, 0, 0))

            def fold_g2e(g2e, fold=8):
                if fold == 1:
                    return g2e
                elif fold == 8:
                    fg2e = np.zeros((npair * (npair + 1) // 2, ), dtype=g2e.dtype)
                    ij = 0
                    ijkl = 0
                    for i in range(self.n_orb):
                        for j in range(0, i + 1):
                            kl = 0
                            for k in range(0, i + 1):
                                for l in range(0, k + 1):
                                    if ij >= kl:
                                        fg2e[ijkl] = g2e[i, j, k, l]
                                        ijkl += 1
                                    kl += 1
                            ij += 1
                    return fg2e
                else:
                    raise ValueError('fold must be 1 or 8')

            fg2e = fold_g2e(self.g2e,fold=8)
            write_g2e(handle, fg2e)
            write_h1e(handle, self.h1e)
            write_epsilon(handle, self.epsilon)
            handle.write(output_format % (self.e_nuc, 0, 0, 0, 0))
    
    def hamiltonian(self, filename):
        """Load fci information from FCIDUMP file.

        Parameters
        ----------
        filename : str
            Name of the fcidump file.
    
        Return
        ------
        """
        from moha.system.operator.operator import Operator,ElectronRepulsion
        from moha.system.operator.hamiltonian import Hamiltonian

        fcidump = self.load(filename)
        n_orb = fcidump.n_orb

        ham = Hamiltonian()
        S = np.identity((n_orb,n_orb))
        Hcore = Operator(fcidump.h1e)
        Eri = ElectronRepulsion(fcidump.g2e)
        ham['S'] = S
        ham['Hcore'] = Hcore
        ham['Eri'] = Eri

        return ham
