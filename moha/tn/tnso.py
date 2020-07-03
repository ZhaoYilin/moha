import numpy as np
import copy
import scipy.linalg
import moha.tn.mpo_consw as mpo_consw
import moha.tn.itools as itools


class TNSO(object):
    """
    """
    def __init__(self,site,N,P={}):
        self.site = site
        self.N=N 
        self.P=P

class MPO(TNSO):
    """
    """
    def __init__(self,site,N,model='Heisenberg',P={'J':1.0,'Jz':1.0,'h':0.5}):
        """
        N: numer of site
        P: dic of all parameters for the model hamilton
           {'J':1.0,'Jz':1.0,'h':1.0}
        """
        self.site = site
        self.N=N 
        self.P=P
        if model=='Ising':
            self.vertices = copy.deepcopy(self.Ising()) 
        elif model=='XXX':
            self.vertices = copy.deepcopy(self.XXX())
        elif model=='XXZ' or model=='Heisenberg':
            self.vertices = copy.deepcopy(self.XXZ())
        elif model=='XYZ':
            self.vertices = copy.deepcopy(self.XYZ())
        elif model=='Hubbard':
            self.vertices = copy.deepcopy(self.Hubbard())
        elif model=='t-J':
            pass

    def Ising(self):
        """
        Ising Hamiltonian
            H = JSzSz - hSx
        P: dic of all parameters for the model hamilton
           {'J':1.0,'h':1.0}
        """
        self.D = 3
        P = self.P
        N = self.N
        site = self.site
        iden = site.operators["id"]
        s_x = site.operators["s_x"]
        s_z = site.operators["s_z"]

        vertices = []

        row = np.zeros((1,2,2,3))
        column = np.zeros((3,2,2,1))
        matrix = np.zeros((3,2,2,3))

        #Left-hand edge is 1x5 matrix
        #[[I,   J*Sz,  -h*Sx]]
        row[0,:, :,0] = iden[:, :]
        row[0,:, :,1] = P['J']*s_z[:, :]
        row[0,:, :,2] = -P['h']*s_x[:, :]

        #Right-hand edge is 5x1 matrix
        #[[-hSx],
        #[Sz],
        #[I]]
        column[0,:, :,0] = -P['h']*s_x[:, :]
        column[1,:, :,0] = s_z[:, :]
        column[2,:, :,0] = iden[:, :]

        #Bulk Hamiltonian MPO
        #[I,    J*Sz,   -h*Sx]
        #[Z,    Z,      Sz]
        #[Z,    Z,      I]
        matrix[0, :, :, :] = row[0, :]
        matrix[:, :, :, 2] = column[:, 0]

        vertices.append(row[:, :])
        for i in range(1, N-1):
            vertices.append(matrix[:, :, :, :])
        vertices.append(column[:, :])

        return vertices

    def XXX(self):
        pass

    def XXZ(self):
        """
        XXZ Hamiltonian
            H = JSzSz - hSx
        Parameters:
            P: parameters of the model Hamiltonians for each term
               J,Jz,h

        """
        print('hello hfwofew')
        self.D =5
        P = self.P
        N = self.N
        site = self.site
        iden = site.operators["id"]
        s_p = site.operators["s_p"]
        s_m = site.operators["s_m"]
        s_z = site.operators["s_z"]

        vertices = []
        
        row = np.zeros((1, 5, 2, 2))
        column = np.zeros((5, 1, 2, 2))
        matrix = np.zeros((5, 5, 2, 2))
        #Left-hand edge is 1x5 matrix
        #[[I,   0.5J*Sp,  0.5J*Sm, Jz*Sz,  h*Sz]]
        row[0, 0, :, :] = iden[:, :]
        row[0, 1, :, :] = P['J'] / 2. * s_p[:, :]
        row[0, 2, :, :] = P['J'] / 2. * s_m[:, :]
        row[0, 3, :, :] = P['Jz'] * s_z[:, :]
        row[0, 4, :, :] = - P['h'] * s_z[:, :]

        #Right-hand edge is 5x1 matrix
        #[[-h*Sz],
        #[Sm],
        #[Sp],
        #[Sz],
        #[I]]
        column[0, 0, :, :] = - P['h'] * s_z[:, :]
        column[1, 0, :, :] = s_m[:, :]
        column[2, 0, :, :] = s_p[:, :]
        column[3, 0, :, :] = s_z[:, :]
        column[4, 0, :, :] = iden[:, :]

        #Bulk Hamiltonian MPO
        #[I,    J*Sp,   J*Sm,   Jz*Sz,  -h*Sz]
        #[Z,    Z,      Z,      Z,      Sz]
        #[Z,    Z,      Z,      Z,      Sm]
        #[Z,    Z,      Z,      Z,      Sp]
        #[Z,    Z,      Z,      Z,      I]
        matrix[0, :, :, :] = row[0, :]
        matrix[:, 4, :, :] = column[:, 0]

        #The whole MPO
        vertices.append(row[:, :])
        for i in range(1, N-1):
            vertices.append(matrix[:, :, :, :])
        vertices.append(column[:, :])

        return vertices


    def XYZ(self):
        pass

    def Hubbard(self):
        """
        XXZ Hamiltonian
            H = JSzSz - hSx
        Parameters:
            P: parameters of the model Hamiltonians for each term
               J,Jz,h

        """
        self.D = 6
        P = self.P
        N = self.N
        site = self.site
        iden = site.operators["id"]
        S_p_up = site.operators["S_p_up"]
        S_p_down = site.operators["S_p_down"]
        S_m_up = site.operators["S_m_up"]
        S_m_down = site.operators["S_m_down"]
        n_up = site.operators["n_up"]
        n_down = site.operators["n_down"]
        F = site.operators["F"]

        vertices = []
        
        row = np.zeros((1, 6, 4, 4))
        column = np.zeros((6, 1, 4, 4))
        matrix = np.zeros((6, 6, 4, 4))
        #Left-hand edge is 1x6 matrix
        #[[I,   S_p_up*F,  S_p_d*F, S_p_d, S_p_up,  U*n_up*n_down]]
        row[0, 0, :, :] = iden[:, :]
        row[0, 1, :, :] = np.dot(S_p_up,F)[:, :]
        row[0, 2, :, :] = np.dot(S_m_up,F)[:, :]
        row[0, 3, :, :] = S_p_down[:, :]
        row[0, 4, :, :] = S_m_down[:, :]
        row[0, 5, :, :] = P['U']*np.dot(n_up,n_down)[:, :]

        #Right-hand edge is 6x1 matrix
        #[[U*n_up*n_down],
        #[-t*S_m_up],
        #[t*S_p_up],
        #[-t*F*S_m_down],
        #[t*F*S_p_down],
        #[I]]
        column[0, 0, :, :] = P['U']*np.dot(n_up,n_down)[:, :]
        column[1, 0, :, :] = -P['t']*S_m_up[:,:]
        column[2, 0, :, :] = P['t']*S_p_up[:,:]
        column[3, 0, :, :] = -P['t']*np.dot(F,S_m_down)[:,:]
        column[4, 0, :, :] = P['t']*np.dot(F,S_p_down)[:,:]
        column[5, 0, :, :] = iden[:, :]

        #Bulk Hamiltonian MPO
        #[I,   S_p_up*F,  S_p_d*F, S_p_d, S_p_up,  U*n_up*n_down]
        #[Z,   Z,         Z,       Z,     Z,       -t*S_m_up]
        #[Z,   Z,         Z,       Z,     Z,       t*S_p_up]
        #[Z,   Z,         Z,       Z,     Z,       -t*F*S_m_down]
        #[Z,   Z,         Z,       Z,     Z,       t*F*S_p_down]
        #[Z,   Z,         Z,       Z,     Z,       I]
        matrix[0, :, :, :] = row[0, :]
        matrix[:, 5, :, :] = column[:, 0]

        #The whole MPO
        vertices.append(row[:, :])
        for i in range(1, N-1):
            vertices.append(matrix[:, :, :, :])
        vertices.append(column[:, :])

        return vertices

    def tJ(self):
        pass

    def one_vertice_exp(self,i,t):
        pass

    def two_vertices_exp(self,i,j,t):
        pass

    def exp(self,t):
        pass


class QCMPO(TNSO):
    """
    """
    def __init__(self,hfham,hfwf,partition,isym=0):
        """
        N: numer of site
        P: dic of all parameters for the model hamilton
           {'J':1.0,'Jz':1.0,'h':1.0}
        """
        h1e,h2e = self.molecular_spin_integral(hfham,hfwf) 
        self.ngroups, self.vertices = copy.deepcopy(self.directHmpo(h1e,h2e,partition,isym=0))
        
    def molecular_spin_integral(self,hfham,hfwf):
        """
        output
            h1e h2e
        """
        C = np.identity(7)
        print(C)
        hfham.operators['electron_repulsion'].basis_transformation(C)
        h1e = hfham.operators['kinetic'].spin
        h1e += hfham.operators['nuclear_attraction'].spin
        h2e = hfham.operators['electron_repulsion'].spin
        return h1e,h2e

    # INPUT: h1e,eri; OUTPUT: MPOs
    def directHmpo(self,h1e,h2e,partition=None,debug=False,iprt=0,isym=0):
        if iprt>0: print('\n[directHmpo]')
        if debug:
            print(' h1e.shape=',h1e.shape)
            print(' h2e.shape=',h2e.shape)
        k = h1e.shape[0]
        if partition == None:
            partition = [[i] for i in range(k)]
        fsupport = self.partFlatten(partition)
        assert len(fsupport) == k
        assert fsupport[-1] == k-1
        ngroups = len(partition)
        # Special case
        if ngroups == 1:
            sl = []
            sc = fsupport
            sr = []
            tmp = mpo_consw.l1r16(h1e,h2e,sl,sc,sr)
            nc1,nc2 = tmp.shape
            tmp = tmp.reshape((1,1,nc1,nc2))
            return ngroups,[tmp]
        ndims = map(lambda x:len(x),partition)
        if debug:
            print(' partition=',partition)
            print(' ngroups  =',ngroups)
            print(' ndims    =',ndims)
            print(' fsupport =',fsupport)
        #
        # LOOP over groups
        #
        tablcr,tabspt,tabdim = self.genTablesLCR(partition,debug)
        wdims = self.genTablesWdims(tabdim,debug)
        wdims = [1]+wdims+[1]
        wfacs = [0]*ngroups
        qnums = [0]*ngroups
        if debug:
            print(' tablcr =',tablcr)
            print(' tabspt =',tabspt)
            print(' tabdim =',tabdim)
            print(' wdims  =',wdims)
        for igroup in range(ngroups):
            diml = wdims[igroup]
            dimr = wdims[igroup+1]
            lg,cg,rg,dims1l,dims1s,dims2l,dims2s = self.genTablesDim12(tabdim,igroup)
            sl,sc,sr = tabspt[igroup]
            wtmp = np.zeros((diml,dimr,2**cg,2**cg))
            a1 = [0]+list(itools.accumulate(dims1l))
            a2 = [0]+list(itools.accumulate(dims2l))
            if debug:
                print('igroup=',igroup,diml,dimr)
                print(' wfacs[i]=',wtmp.shape)
                print(' sl,sc,sr=',sl,sc,sr)
                print(' lg,cg,rg=',(lg,cg,rg))
                print(' dims1l = ',len(dims1l),dims1l)
                print(' dims2l = ',len(dims2l),dims2l)
                print(' accum1 = ',len(a1),a1)
                print(' accum2 = ',len(a2),a2)
            if igroup == 0:
                # row-1: match last index for a1,a2
                wtmp[0:1,a2[0]:a2[1]] = mpo_consw.l1r1(h1e,h2e,sl,sc,sr)
                wtmp[0:1,a2[1]:a2[2]] = mpo_consw.l1r2(h1e,h2e,sl,sc,sr)
                wtmp[0:1,a2[2]:a2[3]] = mpo_consw.l1r3(h1e,h2e,sl,sc,sr)
                if a2[3]<a2[4]:
                    wtmp[0:1,a2[3]:a2[4]] = mpo_consw.l1r4(h1e,h2e,sl,sc,sr)
                if a2[4]<a2[5]:
                    wtmp[0:1,a2[4]:a2[5]] = mpo_consw.l1r5(h1e,h2e,sl,sc,sr)
                wtmp[0:1,a2[6]:a2[7]] = mpo_consw.l1r7(h1e,h2e,sl,sc,sr)
                wtmp[0:1,a2[10]:a2[11]] = mpo_consw.l1r11(h1e,h2e,sl,sc,sr)
                wtmp[0:1,a2[12]:a2[13]] = mpo_consw.l1r13(h1e,h2e,sl,sc,sr)
                wtmp[0:1,a2[14]:a2[15]] = mpo_consw.l1r15(h1e,h2e,sl,sc,sr)
                wtmp[0:1,a2[15]:a2[16]] = mpo_consw.l1r16(h1e,h2e,sl,sc,sr)
                #----------------------------------------------
                # Reordering of Qterm
                ijdx = self.sortQlc(lg,cg)+a2[7]
                assert len(ijdx)==a2[11]-a2[7]
                wtmp[:,a2[7]:a2[11]] = wtmp[:,ijdx].copy()
                #----------------------------------------------
            elif igroup == ngroups-1:
                # col-1
                wtmp[a1[0]:a1[1]  ,0:1] = mpo_consw.l1r16(h1e,h2e,sl,sc,sr)
                wtmp[a1[1]:a1[2]  ,0:1] = mpo_consw.l2r16(h1e,h2e,sl,sc,sr)
                wtmp[a1[3]:a1[4]  ,0:1] = mpo_consw.l4r16(h1e,h2e,sl,sc,sr)
                if a1[5] < a1[6]:
                    wtmp[a1[5]:a1[6],0:1] = mpo_consw.l6r16(h1e,h2e,sl,sc,sr)
                if a1[8] < a1[9]:
                    wtmp[a1[8]:a1[9],0:1] = mpo_consw.l9r16(h1e,h2e,sl,sc,sr)
                wtmp[a1[11]:a1[12],0:1] = mpo_consw.l12r16(h1e,h2e,sl,sc,sr)
                wtmp[a1[12]:a1[13],0:1] = mpo_consw.l13r16(h1e,h2e,sl,sc,sr)
                wtmp[a1[13]:a1[14],0:1] = mpo_consw.l14r16(h1e,h2e,sl,sc,sr)
                wtmp[a1[14]:a1[15],0:1] = mpo_consw.l15r16(h1e,h2e,sl,sc,sr)
                wtmp[a1[15]:a1[16],0:1] = mpo_consw.l16r16(h1e,h2e,sl,sc,sr)
                #----------------------------------------------
                # Reordering of Pterm
                ijdx = self.sortPcr(cg,rg)
                assert len(ijdx)==a1[8]-a1[5]
                wtmp[a1[5]:a1[8],:] = wtmp[ijdx+a1[5],:].copy()
                assert len(ijdx)==a1[11]-a1[8]
                wtmp[a1[8]:a1[11],:] = wtmp[ijdx+a1[8],:].copy()
                #----------------------------------------------
            else:
                # row-1 => H
                wtmp[a1[0]:a1[1],a2[0]:a2[1]] = mpo_consw.l1r1(h1e,h2e,sl,sc,sr)
                wtmp[a1[0]:a1[1],a2[1]:a2[2]] = mpo_consw.l1r2(h1e,h2e,sl,sc,sr)
                wtmp[a1[0]:a1[1],a2[2]:a2[3]] = mpo_consw.l1r3(h1e,h2e,sl,sc,sr)
                if a2[3]<a2[4]:
                    wtmp[a1[0]:a1[1],a2[3]:a2[4]] = mpo_consw.l1r4(h1e,h2e,sl,sc,sr)
                if a2[4]<a2[5]:
                    wtmp[a1[0]:a1[1],a2[4]:a2[5]] = mpo_consw.l1r5(h1e,h2e,sl,sc,sr)
                wtmp[a1[0]:a1[1],a2[6]:a2[7]] = mpo_consw.l1r7(h1e,h2e,sl,sc,sr)
                wtmp[a1[0]:a1[1],a2[10]:a2[11]] = mpo_consw.l1r11(h1e,h2e,sl,sc,sr)
                wtmp[a1[0]:a1[1],a2[12]:a2[13]] = mpo_consw.l1r13(h1e,h2e,sl,sc,sr)
                wtmp[a1[0]:a1[1],a2[14]:a2[15]] = mpo_consw.l1r15(h1e,h2e,sl,sc,sr)
                wtmp[a1[0]:a1[1],a2[15]:a2[16]] = mpo_consw.l1r16(h1e,h2e,sl,sc,sr)
                # row-2 => a+[i]
                wtmp[a1[1]:a1[2],a2[15]:a2[16]] = mpo_consw.l2r16(h1e,h2e,sl,sc,sr)
                # row-3 => a+[i+1,K]
                wtmp[a1[2]:a1[3],a2[1]:a2[2]] = mpo_consw.l3r2(h1e,h2e,sl,sc,sr)
                # row-4 => a[i]
                wtmp[a1[3]:a1[4],a2[15]:a2[16]] = mpo_consw.l4r16(h1e,h2e,sl,sc,sr)
                # row-5 => a[i+1,K]
                wtmp[a1[4]:a1[5],a2[2]:a2[3]] = mpo_consw.l5r3(h1e,h2e,sl,sc,sr)
                # row-6 => A+[i][i]
                if a1[5] < a1[6]:
                    wtmp[a1[5]:a1[6],a2[15]:a2[16]] = mpo_consw.l6r16(h1e,h2e,sl,sc,sr)
                # row-7 => A+[i][i+1,K]
                wtmp[a1[6]:a1[7],a2[1]:a2[2]] = mpo_consw.l7r2(h1e,h2e,sl,sc,sr)
                # row-8 => A+[i+1,K][i+1,K]
                if a1[7] < a1[8]:
                    wtmp[a1[7]:a1[8],a2[3]:a2[4]] = mpo_consw.l8r4(h1e,h2e,sl,sc,sr)
                # row-9 => A[i][i]
                if a1[8] < a1[9]:
                    wtmp[a1[8]:a1[9],a2[15]:a2[16]] = mpo_consw.l9r16(h1e,h2e,sl,sc,sr)
                # row-10 => A[i][i+1,K]
                wtmp[a1[9]:a1[10],a2[2]:a2[3]] = mpo_consw.l10r3(h1e,h2e,sl,sc,sr)
                # row-11 => A[i+1,K]A[i+1,K]
                if a1[10]<a1[11]:
                    wtmp[a1[10]:a1[11],a2[4]:a2[5]] = mpo_consw.l11r5(h1e,h2e,sl,sc,sr)
                # row-12 => S[1,i-1]
                wtmp[a1[11]:a1[12],a2[5]:a2[6]] = mpo_consw.l12r6(h1e,h2e,sl,sc,sr)
                wtmp[a1[11]:a1[12],a2[15]:a2[16]] = mpo_consw.l12r16(h1e,h2e,sl,sc,sr)
                # row-13 => Q[1,i-1][1,i-1]
                wtmp[a1[12]:a1[13],a2[1]:a2[2]] = mpo_consw.l13r2(h1e,h2e,sl,sc,sr)
                wtmp[a1[12]:a1[13],a2[2]:a2[3]] = mpo_consw.l13r3(h1e,h2e,sl,sc,sr)
                wtmp[a1[12]:a1[13],a2[7]:a2[8]] = mpo_consw.l13r8(h1e,h2e,sl,sc,sr)
                wtmp[a1[12]:a1[13],a2[15]:a2[16]] = mpo_consw.l13r16(h1e,h2e,sl,sc,sr)
                # row-14 => T1[1,i-1]
                wtmp[a1[13]:a1[14],a2[1]:a2[2]] = mpo_consw.l14r2(h1e,h2e,sl,sc,sr)
                wtmp[a1[13]:a1[14],a2[2]:a2[3]] = mpo_consw.l14r3(h1e,h2e,sl,sc,sr)
                if a2[4]<a2[5]:
                    wtmp[a1[13]:a1[14],a2[4]:a2[5]] = mpo_consw.l14r5(h1e,h2e,sl,sc,sr)
                wtmp[a1[13]:a1[14],a2[8]:a2[9]] = mpo_consw.l14r9(h1e,h2e,sl,sc,sr)
                wtmp[a1[13]:a1[14],a2[11]:a2[12]] = mpo_consw.l14r12(h1e,h2e,sl,sc,sr)
                wtmp[a1[13]:a1[14],a2[15]:a2[16]] = mpo_consw.l14r16(h1e,h2e,sl,sc,sr)
                # row-15 => T3[1,i-1]
                wtmp[a1[14]:a1[15],a2[1]:a2[2]] = mpo_consw.l15r2(h1e,h2e,sl,sc,sr)
                wtmp[a1[14]:a1[15],a2[2]:a2[3]] = mpo_consw.l15r3(h1e,h2e,sl,sc,sr)
                if a2[3]<a2[4]:
                    wtmp[a1[14]:a1[15],a2[3]:a2[4]] = mpo_consw.l15r4(h1e,h2e,sl,sc,sr)
                wtmp[a1[14]:a1[15],a2[9]:a2[10]] = mpo_consw.l15r10(h1e,h2e,sl,sc,sr)
                wtmp[a1[14]:a1[15],a2[13]:a2[14]] = mpo_consw.l15r14(h1e,h2e,sl,sc,sr)
                wtmp[a1[14]:a1[15],a2[15]:a2[16]] = mpo_consw.l15r16(h1e,h2e,sl,sc,sr)
                # row-16 => I[i,K]
                wtmp[a1[15]:a1[16],a2[15]:a2[16]] = mpo_consw.l16r16(h1e,h2e,sl,sc,sr)
                #----------------------------------------------
                # Reordering of Qterm
                ijdx = self.sortQlc(lg,cg)+a2[7]
                assert len(ijdx)==a2[11]-a2[7]
                wtmp[:,a2[7]:a2[11]] = wtmp[:,ijdx].copy()
                #----------------------------------------------
                # Reordering of Pterm
                ijdx = self.sortPcr(cg,rg)
                assert len(ijdx)==a1[8]-a1[5]
                wtmp[a1[5]:a1[8],:] = wtmp[ijdx+a1[5],:].copy()
                assert len(ijdx)==a1[11]-a1[8]
                wtmp[a1[8]:a1[11],:] = wtmp[ijdx+a1[8],:].copy()
            #----------------------------------------------
            # Store
            wfacs[igroup] = wtmp.copy()
            #----------------------------------------------
            # Qnumbers
            #----------------------------------------------
            qtmp = [None]*dimr
            if isym == 1:
                qtmp[a2[0] :a2[1]]  = [ [0]  ]*(a2[1] -a2[0] )
                qtmp[a2[1] :a2[2]]  = [ [-1] ]*(a2[2] -a2[1] )
                qtmp[a2[2] :a2[3]]  = [ [+1] ]*(a2[3] -a2[2] )
                qtmp[a2[3] :a2[4]]  = [ [-2] ]*(a2[4] -a2[3] )
                qtmp[a2[4] :a2[5]]  = [ [+2] ]*(a2[5] -a2[4] )
                qtmp[a2[5] :a2[7]]  = [ [-1] ]*(a2[7] -a2[5] )
                qtmp[a2[7] :a2[11]] = [ [0]  ]*(a2[11]-a2[7] )
                qtmp[a2[11]:a2[13]] = [ [+1] ]*(a2[13]-a2[11])
                qtmp[a2[13]:a2[15]] = [ [-1] ]*(a2[15]-a2[13])
                qtmp[a2[15]:a2[16]] = [ [0]  ]*(a2[16]-a2[15])
                if igroup == ngroups-1:
                    qnums[igroup] = [[0]]
                else:
                    qnums[igroup] = copy.deepcopy(qtmp)
            elif isym == 2:
                # Only works for spin-orbital case
                assert k == ngroups
                if igroup%2 == 0:
                    sz = 0.5
                else:
                    sz = -0.5
                qtmp[a2[0] :a2[1]]  = [ [0,0]     ]*(a2[1] -a2[0] )
                qtmp[a2[1] :a2[2]]  = [ [-1,-sz]  ]*(a2[2] -a2[1] ) # T2
                qtmp[a2[2] :a2[3]]  = [ [+1,sz]   ]*(a2[3] -a2[2] ) # S+T4
                qtmp[a2[3] :a2[4]]  = [ [-2,-2*sz]]*(a2[4] -a2[3] ) # P2ann: In fact, they will
                qtmp[a2[4] :a2[5]]  = [ [+2,2*sz] ]*(a2[5] -a2[4] ) # P2cre: not exisit at all.
                qtmp[a2[5] :a2[7]]  = [ [-1,-sz]  ]*(a2[7] -a2[5] )
                qtmp[a2[7] :a2[11]] = [ [0,0]     ]*(a2[11]-a2[7] )
                qtmp[a2[11]:a2[13]] = [ [+1,sz]   ]*(a2[13]-a2[11])
                qtmp[a2[13]:a2[15]] = [ [-1,-sz]  ]*(a2[15]-a2[13])
                qtmp[a2[15]:a2[16]] = [ [0,0]     ]*(a2[16]-a2[15])
                if igroup == ngroups-1:
                    qnums[igroup] = [[0,0]]
                else:
                    qnums[igroup] = copy.deepcopy(qtmp)
        #----------------------------------------------
        # Finally, form MPO
        if iprt>0: print(' wdims=',wdims)
        if isym == 0:
            return ngroups,wfacs
        else:
            return ngroups,wfacs,qnums

    def antisymmetrizeTwoBodyOpers(self,eri):
        sbas = eri.shape[0]
        eriA=eri-eri.transpose(0,3,2,1)
        eriA=eriA.transpose(0,2,1,3).copy()
        # Unique representation: V[ijkl] (i<j,k<l) = -<ij||kl>
        aeri=np.zeros(eriA.shape)
        for j in range(sbas):
            for i in range(j):
                for l in range(sbas):
                    for k in range(l):
                        aeri[i,j,k,l] = -eriA[i,j,k,l]
        return aeri

    def partFlatten(self,partition):
        support = []
        for ipart in partition:
            support += ipart
        return support

    def genTablesLCR(self,partition,debug=False):
        ngroups = len(partition)
        tablcr = []
        tabspt = []
        tabdim = []
        for igroup in range(ngroups):
            # Partition of the basis
            lg = partition[:igroup]
            cg = [partition[igroup]]
            rg = partition[igroup+1:]
            sl = self.partFlatten(lg)
            sc = self.partFlatten(cg)
            sr = self.partFlatten(rg)
            # cardinality of basis
            nl = len(sl)
            nc = len(sc)
            nr = len(sr)
            tablcr.append([lg,cg,rg])
            tabspt.append([sl,sc,sr])
            tabdim.append([nl,nc,nr])
            if debug:
                print(' igroup =',igroup)
                print(' lg = ',lg)
                print(' cg = ',cg)
                print(' rg = ',rg)
                print(' sl = ',sl)
                print(' sc = ',sc)
                print(' sr = ',sr)
                print(' nl,nc,nr = ',nl,nc,nr)
        return tablcr,tabspt,tabdim

    def genTablesDim12(self,tabdim,igroup):
        lg,cg,rg = tabdim[igroup]
        lg2 = int(lg*(lg-1)/2)
        cg2 = int(cg*(cg-1)/2)
        rg2 = int(rg*(rg-1)/2)
        dims1l = [1,cg,rg,cg,rg,cg2,cg*rg,rg2,cg2,cg*rg,rg2,lg,lg**2,lg,lg,1]
        dims1s = [1,cg+rg,cg+rg,cg2+cg*rg+rg2,cg2+cg*rg+rg2,lg,lg**2,lg,lg,1]
        dims2l = [1,rg,rg,rg2,rg2,lg,cg,lg**2,lg*cg,lg*cg,cg**2,lg,cg,lg,cg,1]
        dims2s = [1,rg,rg,rg2,rg2,lg+cg,lg**2+lg*cg+lg*cg+cg**2,lg+cg,lg+cg,1]
        return lg,cg,rg,dims1l,dims1s,dims2l,dims2s

    def genTablesWdims(self,tabdim,debug=False):
        #
        # W[a1]*W[a1]
        # W[a1]*W[a1,a2]*W[a2]
        # W[a1]*W[a1,a2]*W[a2,a3]*W[a3]
        #
        ngroups = len(tabdim)
        wdims = []
        for igroup in range(ngroups):
            lg,cg,rg,dims1l,dims1s,dims2l,dims2s = self.genTablesDim12(tabdim,igroup)
            dim1 = np.sum(dims1l)
            dim2 = np.sum(dims2l)
            wdims.append(dim2)
            if debug:
                print('igroup',igroup,'dimlcr=',(lg,cg,rg))
                print('dim[%d]='%(igroup-1),dim1,'dim[%d]='%igroup,dim2)
                print(' d1l=',dims1l)
                print(' d1s=',dims1s)
                print(' d2l=',dims2l)
                print(' d2s=',dims2s)
        wdims.pop()
        assert len(wdims) == ngroups-1
        if debug: print('Final wdims=',wdims)
        return wdims

    def sortPcr(self,cg,rg):
        def id(i,j,n):
            return i*(2*n-i-3)/2+j-1
        n = cg+rg
        ijcc = [id(i,j,n) for i in range(cg) for j in range(cg) if i<j]
        ijcr = [id(i,j+cg,n) for i in range(cg) for j in range(rg)]
        ijrr = [id(i+cg,j+cg,n) for i in range(rg) for j in range(rg) if i<j]
        ij = ijcc+ijcr+ijrr
        ij = np.array(ij)
        ijdx = np.argsort(ij)
        return ijdx

    def sortQlc(self,lg,cg):
        def id(i,j,n):
            return i*n+j
        n = lg+cg
        ijll = [id(i,j,n) for i in range(lg) for j in range(lg)]
        ijlc = [id(i,j+lg,n) for i in range(lg) for j in range(cg)]
        ijcl = [id(i+lg,j,n) for i in range(cg) for j in range(lg)]
        ijcc = [id(i+lg,j+lg,n) for i in range(cg) for j in range(cg)]
        ij = ijll+ijlc+ijcl+ijcc
        ij = np.array(ij)
        ijdx = np.argsort(ij)
        return ijdx

