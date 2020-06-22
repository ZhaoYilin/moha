#
# Direct construciton of MPO for H or T2 via W factors
#
import copy
import numpy
import scipy.linalg
import mpo_consw
import itools

# INPUT: h1e,eri; OUTPUT: MPOs
def directHmpo(h1e,h2e,partition=None,debug=False,iprt=0,isym=0):
    if iprt>0: print('\n[directHmpo]')
    if debug:
        print(' h1e.shape=',h1e.shape)
        print(' h2e.shape=',h2e.shape)
    k = h1e.shape[0]
    if partition == None:
        partition = [[i] for i in range(k)]
    fsupport = partFlatten(partition)
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
    tablcr,tabspt,tabdim = genTablesLCR(partition,debug)
    wdims = genTablesWdims(tabdim,debug)
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
        lg,cg,rg,dims1l,dims1s,dims2l,dims2s = genTablesDim12(tabdim,igroup)
        sl,sc,sr = tabspt[igroup]
        wtmp = numpy.zeros((diml,dimr,2**cg,2**cg))
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
            ijdx = sortQlc(lg,cg)+a2[7]
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
            ijdx = sortPcr(cg,rg)
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
            ijdx = sortQlc(lg,cg)+a2[7]
            assert len(ijdx)==a2[11]-a2[7]
            wtmp[:,a2[7]:a2[11]] = wtmp[:,ijdx].copy()
            #----------------------------------------------
            # Reordering of Pterm
            ijdx = sortPcr(cg,rg)
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

def antisymmetrizeTwoBodyOpers(eri):
    sbas = eri.shape[0]
    eriA=eri-eri.transpose(0,3,2,1)
    eriA=eriA.transpose(0,2,1,3).copy()
    # Unique representation: V[ijkl] (i<j,k<l) = -<ij||kl>
    aeri=numpy.zeros(eriA.shape)
    for j in range(sbas):
        for i in range(j):
            for l in range(sbas):
                for k in range(l):
                    aeri[i,j,k,l] = -eriA[i,j,k,l]
    return aeri

def partFlatten(partition):
    support = []
    for ipart in partition:
        support += ipart
    return support

def genTablesLCR(partition,debug=False):
    ngroups = len(partition)
    tablcr = []
    tabspt = []
    tabdim = []
    for igroup in range(ngroups):
        # Partition of the basis
        lg = partition[:igroup]
        cg = [partition[igroup]]
        rg = partition[igroup+1:]
        sl = partFlatten(lg)
        sc = partFlatten(cg)
        sr = partFlatten(rg)
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

def genTablesDim12(tabdim,igroup):
    lg,cg,rg = tabdim[igroup]
    lg2 = int(lg*(lg-1)/2)
    cg2 = int(cg*(cg-1)/2)
    rg2 = int(rg*(rg-1)/2)
    dims1l = [1,cg,rg,cg,rg,cg2,cg*rg,rg2,cg2,cg*rg,rg2,lg,lg**2,lg,lg,1]
    dims1s = [1,cg+rg,cg+rg,cg2+cg*rg+rg2,cg2+cg*rg+rg2,lg,lg**2,lg,lg,1]
    dims2l = [1,rg,rg,rg2,rg2,lg,cg,lg**2,lg*cg,lg*cg,cg**2,lg,cg,lg,cg,1]
    dims2s = [1,rg,rg,rg2,rg2,lg+cg,lg**2+lg*cg+lg*cg+cg**2,lg+cg,lg+cg,1]
    return lg,cg,rg,dims1l,dims1s,dims2l,dims2s

def genTablesWdims(tabdim,debug=False):
    #
    # W[a1]*W[a1]
    # W[a1]*W[a1,a2]*W[a2]
    # W[a1]*W[a1,a2]*W[a2,a3]*W[a3]
    #
    ngroups = len(tabdim)
    wdims = []
    for igroup in range(ngroups):
        lg,cg,rg,dims1l,dims1s,dims2l,dims2s = genTablesDim12(tabdim,igroup)
        dim1 = numpy.sum(dims1l)
        dim2 = numpy.sum(dims2l)
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

def sortPcr(cg,rg):
    def id(i,j,n):
        return i*(2*n-i-3)/2+j-1
    n = cg+rg
    ijcc = [id(i,j,n) for i in range(cg) for j in range(cg) if i<j]
    ijcr = [id(i,j+cg,n) for i in range(cg) for j in range(rg)]
    ijrr = [id(i+cg,j+cg,n) for i in range(rg) for j in range(rg) if i<j]
    ij = ijcc+ijcr+ijrr
    ij = numpy.array(ij)
    ijdx = numpy.argsort(ij)
    return ijdx

def sortQlc(lg,cg):
    def id(i,j,n):
        return i*n+j
    n = lg+cg
    ijll = [id(i,j,n) for i in range(lg) for j in range(lg)]
    ijlc = [id(i,j+lg,n) for i in range(lg) for j in range(cg)]
    ijcl = [id(i+lg,j,n) for i in range(cg) for j in range(lg)]
    ijcc = [id(i+lg,j+lg,n) for i in range(cg) for j in range(cg)]
    ij = ijll+ijlc+ijcl+ijcc
    ij = numpy.array(ij)
    ijdx = numpy.argsort(ij)
    return ijdx

#==================
# Test subroutines
#==================
if __name__ == '__main__':

    def genRandomH(k):
        numpy.random.seed(1)
        h1e = numpy.random.uniform(-1,1,size=(k,k))
        h2e = numpy.random.uniform(-1,1,size=(k,k,k,k))
        h1e = 0.5*(h1e+h1e.T)
        h2e = 0.5*(h2e+h2e.transpose((2,3,0,1)))
        h2e = antisymmetrizeTwoBodyOpers(h2e)
        return h1e,h2e

    k = 10
    h1e,h2e = genRandomH(k)
    partition = [[0,1,2,3],[4,5],[6,7],[8,9]]
    ngroups,wfacs = directHmpo(h1e,h2e,partition,isym=0,debug=True)
    print(ngroups)
    for i in range(ngroups):
        print(wfacs[i].shape)
