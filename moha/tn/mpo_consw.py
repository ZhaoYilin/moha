#
# Construct Wfactors for a given Hamiltonian (H1,H2)
#
import numpy
import moha.tn.jwtrans as jwtrans
from functools import reduce


#======
# AUX
#======
def idnC(sc):
    nc = len(sc)
    if nc == 0:
        print('error: mpo_consw.idnC')
        exit()
    else:
        idn = jwtrans.idn
        for i in range(nc-1):
            idn = numpy.kron(idn,jwtrans.idn)
    return idn

def sgnC(sc):
    nc = len(sc)
    if nc == 0:
        print('error: mpo_consw.sgnC')
        exit()
    else:
        sgn = jwtrans.sgn
        for i in range(nc-1):
            sgn = numpy.kron(sgn,jwtrans.sgn)
    return sgn

def sqC(sc,ic,iop):
    nc = len(sc)
    if nc == 0:
        print('error: mpo_consw.sqC')
        exit()
    else:
        ops = [0]*nc
        if iop == 1: 
            ops[ic] = jwtrans.cre
        elif iop == 0:
            ops[ic] = jwtrans.ann
        for i in range(ic):
            ops[i] = jwtrans.sgn
        for i in range(ic+1,nc):
            ops[i] = jwtrans.idn
        tmp = reduce(numpy.kron,ops)
    return tmp

def eint(a,b,c,d):
    if a > b: ab = a*(a+1)/2 + b
    else: ab = b*(b+1)/2 + a
    if c > d: cd = c*(c+1)/2 + d
    else: cd = d*(d+1)/2 + c
    if ab > cd: abcd = ab*(ab+1)/2 + cd
    else: abcd = cd*(cd+1)/2 + ab
    return int(abcd)
#======
# Main
#======
# dims1l = [1,cg,rg,cg,rg,cg2,cg*rg,rg2,cg2,cg*rg,rg2,lg,lg**2,lg,lg,1]
# dims2l = [1,rg,rg,rg2,rg2,lg,cg,lg**2,lg*cg,lg*cg,cg**2,lg,cg,lg,cg,1]
#======

#------------
# row-1 => H
#------------
# w[1,1]=I
def l1r1(h1e,h2e,sl,sc,sr):
    tmp = idnC(sc)
    return tmp

# w[1,2]=T2
def l1r2(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nr,2**nc,2**nc))
    for idxp,orbp in enumerate(sc):
        op_p = sqC(sc,idxp,1)
        for idxr,orbr in enumerate(sc):
            op_r = sqC(sc,idxr,0)
            for idxs,orbs in enumerate(sc):
                op_s = sqC(sc,idxs,0)
                # r < s
                if orbr < orbs:
                    op = op_p.dot(op_r.dot(op_s))
                    for idxq,orbq in enumerate(sr):
                        tmp[idxq] += h2e[eint(orbp,orbq,orbr,orbs)]*op
    sgn = sgnC(sc)
    for idxq,orbq in enumerate(sr):
        tmp[idxq] = tmp[idxq].dot(sgn)
    return tmp

# w[1,3]=Sq+T4s
def l1r3(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nr,2**nc,2**nc))
    for idxp,orbp in enumerate(sc):
        op_p = sqC(sc,idxp,1)
        # +h1e
        for idxq,orbq in enumerate(sr):
            tmp[idxq] += h1e[orbp,orbq]*op_p
        for idxq,orbq in enumerate(sc):
            op_q = sqC(sc,idxq,1)
            if orbp < orbq:
                for idxr,orbr in enumerate(sc):
                    op_r = sqC(sc,idxr,0)
                    op = op_p.dot(op_q.dot(op_r))
                    for idxs,orbs in enumerate(sr):
                        tmp[idxs] += h2e[eint(orbp,orbq,orbr,orbs)]*op
    sgn = sgnC(sc)
    for idxq,orbq in enumerate(sr):
        tmp[idxq] = tmp[idxq].dot(sgn)
    return tmp

# w[1,4]=P
def l1r4(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    if int(nr*(nr-1)/2)==0:
        print('error: mpo_consw.l1r4')
        exit()
    tmp = numpy.zeros((int(nr*(nr-1)/2),2**nc,2**nc))
    for idxr,orbr in enumerate(sc):
        op_r = sqC(sc,idxr,0)
        for idxs,orbs in enumerate(sc):
            op_s = sqC(sc,idxs,0)
            if orbr < orbs:
                op = op_r.dot(op_s)
                ipq = 0
                for idxp,orbp in enumerate(sr):
                    for idxq,orbq in enumerate(sr):
                        if orbp < orbq:
                            tmp[ipq] += h2e[eint(orbp,orbq,orbr,orbs)]*op
                            ipq += 1
    return tmp

# w[1,5]=P^+
def l1r5(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    if int(nr*(nr-1)/2)==0:
        print('error: mpo_consw.l1r5')
        exit()
    tmp = numpy.zeros((int(nr*(nr-1)/2),2**nc,2**nc))
    for idxr,orbr in enumerate(sc):
        op_r = sqC(sc,idxr,1)
        for idxs,orbs in enumerate(sc):
            op_s = sqC(sc,idxs,1)
            if orbr < orbs:
                op = op_r.dot(op_s)
                ipq = 0
                for idxp,orbp in enumerate(sr):
                    for idxq,orbq in enumerate(sr):
                        if orbp < orbq:
                            tmp[ipq] += h2e[eint(orbr,orbs,orbp,orbq)]*op
                            ipq += 1
    return tmp

# w[1,7]=-tilde{a}
def l1r7(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nc,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxr,orbr in enumerate(sc):
        tmp[idxr] = -numpy.dot(sqC(sc,idxr,0),sgn)
    return tmp

# w[1,11]=(-a^+a)
def l1r11(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nc*nc,2**nc,2**nc))
    irs = 0
    for idxr,orbr in enumerate(sc):
        for idxs,orbs in enumerate(sc):
            tmp[irs] = -numpy.dot(sqC(sc,idxr,1),sqC(sc,idxs,0))
            irs +=1
    return tmp

# w[1,13]=tilde{a}^+
def l1r13(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nc,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxr,orbr in enumerate(sc):
        tmp[idxr] = numpy.dot(sqC(sc,idxr,1),sgn)
    return tmp

# w[1,15]=tilde{a}
def l1r15(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nc,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxr,orbr in enumerate(sc):
        tmp[idxr] = numpy.dot(sqC(sc,idxr,0),sgn)
    return tmp

# w[1,16]=Hc
def l1r16(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((2**nc,2**nc))
    # p^+q
    for idxp,orbp in enumerate(sc):
        op_p = sqC(sc,idxp,1)
        for idxq,orbq in enumerate(sc):
            op_q = sqC(sc,idxq,0)
            op = op_p.dot(op_q)	 
            tmp += h1e[orbp,orbq]*op
    # p^+q^+rs
    for idxp,orbp in enumerate(sc):
        op_p = sqC(sc,idxp,1)
        for idxq,orbq in enumerate(sc):
            op_q = sqC(sc,idxq,1)
            if orbp < orbq:
                op_pq = op_p.dot(op_q)	 
                for idxr,orbr in enumerate(sc):
                    op_r = sqC(sc,idxr,0)
                    for idxs,orbs in enumerate(sc):
                        op_s = sqC(sc,idxs,0)
                        if orbr < orbs:
                            op_rs = op_r.dot(op_s)		  
                            op = op_pq.dot(op_rs)
                            tmp += h2e[eint(orbp,orbq,orbr,orbs)]*op
    return tmp

#----------------
# row-2 => a+[i]
#----------------
# w[2,16]=a+
def l2r16(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nc,1,2**nc,2**nc))
    for idxr,orbr in enumerate(sc):
        tmp[idxr,0] = sqC(sc,idxr,1)
    return tmp

#--------------------
# row-3 => a+[i+1,K]
#--------------------
# w[3,2]=tI
def l3r2(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nr,nr,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxr,orbr in enumerate(sr):
        tmp[idxr,idxr] = sgn.copy()
    return tmp

#---------------
# row-4 => a[i]
#---------------
# w[4,16]=a
def l4r16(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nc,1,2**nc,2**nc))
    for idxr,orbr in enumerate(sc):
        tmp[idxr,0] = sqC(sc,idxr,0)
    return tmp

#--------------------
# row-5 => a[i+1,K]
#--------------------
# w[5,3]=tI
def l5r3(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nr,nr,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxr,orbr in enumerate(sr):
        tmp[idxr,idxr] = sgn.copy()
    return tmp

#--------------------
# row-6 => A+[i][i]
#--------------------
# w[6,16]=(a^+a^+)
def l6r16(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    if int(nc*(nc-1)/2)==0:
        print('error: mpo_consw.l6r16')
        exit()
    tmp = numpy.zeros((int(nc*(nc-1)/2),1,2**nc,2**nc))
    irs = 0
    for idxr,orbr in enumerate(sc):
        for idxs,orbs in enumerate(sc):
            if orbr < orbs:
                tmp[irs,0] = numpy.dot(sqC(sc,idxr,1),sqC(sc,idxs,1))
                irs +=1
    return tmp

#-----------------------
# row-7 => A+[i][i+1,K]
#-----------------------
# w[7,2]=(ta^+)
def l7r2(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nc,nr,nr,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxr,orbr in enumerate(sc):
        for idxs,orbs in enumerate(sr):
            tmp[idxr,idxs,idxs] = numpy.dot(sqC(sc,idxr,1),sgn)
    tmp = tmp.reshape((nc*nr,nr,2**nc,2**nc))
    return tmp

#---------------------------
# row-8 => A+[i+1,K][i+1,K]
#---------------------------
# w[8,4]=I
def l8r4(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    nr2 = int(nr*(nr-1)/2)
    if nr2==0:
        print('error: mpo_consw.l8r4')
        exit()
    tmp = numpy.zeros((nr2,nr2,2**nc,2**nc))
    idn = idnC(sc)
    for irs in range(nr2):
        tmp[irs,irs] = idn.copy()
    return tmp

#------------------
# row-9 => A[i][i]
#------------------
# w[9,16]=(a*a)
def l9r16(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    if int(nc*(nc-1)/2)==0:
        print('error: mpo_consw.l9r16')
        exit()
    tmp = numpy.zeros((int(nc*(nc-1)/2),1,2**nc,2**nc))
    irs = 0
    for idxr,orbr in enumerate(sc):
        for idxs,orbs in enumerate(sc):
            if orbr < orbs:
                tmp[irs,0] = numpy.dot(sqC(sc,idxr,0),sqC(sc,idxs,0))
                irs +=1
    return tmp

#-----------------------
# row-10 => A[i][i+1,K]
#-----------------------
# w[10,3]=ta
def l10r3(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nc,nr,nr,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxr,orbr in enumerate(sc):
        for idxs,orbs in enumerate(sr):
            tmp[idxr,idxs,idxs] = numpy.dot(sqC(sc,idxr,0),sgn)
    tmp = tmp.reshape((nc*nr,nr,2**nc,2**nc))
    return tmp

#----------------------------
# row-11 => A[i+1,K]A[i+1,K]
#----------------------------
# w[11,5]=I
def l11r5(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    nr2 = int(nr*(nr-1)/2)
    if nr2==0:
        print('error: mpo_consw.l11r5')
        exit()
    tmp = numpy.zeros((nr2,nr2,2**nc,2**nc))
    idn = idnC(sc)
    for irs in range(nr2):
        tmp[irs,irs] = idn.copy()
    return tmp

#--------------------
# row-12 => S[1,i-1]
#--------------------
# w[12,6]=tI
def l12r6(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nl,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxr,orbr in enumerate(sl):
        tmp[idxr,idxr] = sgn.copy()
    return tmp

# w[12,16]=S
def l12r16(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,1,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxp,orbp in enumerate(sc):
        op_p = sqC(sc,idxp,1)
        for idxq,orbq in enumerate(sl):
            tmp[idxq,0] += h1e[orbp,orbq]*op_p
    return tmp

#---------------------------
# row-13 => Q[1,i-1][1,i-1]
#---------------------------
# w[13,2]=-v[pL,qR,rL,sC]ta_sC
def l13r2(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nl,nr,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxs,orbs in enumerate(sc):
        op_s = sqC(sc,idxs,0)
        op = op_s.dot(sgn)
        for idxp,orbp in enumerate(sl):
            for idxr,orbr in enumerate(sl):
                for idxq,orbq in enumerate(sr):
                    tmp[idxp,idxr,idxq] += -h2e[eint(orbp,orbq,orbr,orbs)]*op
    tmp = tmp.reshape((nl*nl,nr,2**nc,2**nc))
    return tmp

# w[13,3]=v[pL,qC,rL,sR]*ta^+_qC
def l13r3(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nl,nr,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxq,orbq in enumerate(sc):
        op_q = sqC(sc,idxq,1)
        op = op_q.dot(sgn)
        for idxp,orbp in enumerate(sl):
            for idxr,orbr in enumerate(sl):
                for idxs,orbs in enumerate(sr):
                    tmp[idxp,idxr,idxs] += h2e[eint(orbp,orbq,orbr,orbs)]*op
    tmp = tmp.reshape((nl*nl,nr,2**nc,2**nc))
    return tmp

# w[13,8]=I
#-------------------------------
# It should be understood as
# Q[pq] = M[pq,rs]*Q[rs]
# where M[pq,rs]=delta(pq,rs)*I
#-------------------------------
def l13r8(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl*nl,nl*nl,2**nc,2**nc))
    idn = idnC(sc)
    for ipr in range(nl*nl):
        tmp[ipr,ipr] = idn.copy()
    return tmp

# w[13,16]=Q
def l13r16(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nl,1,2**nc,2**nc))
    for idxq,orbq in enumerate(sc):
        op_q = sqC(sc,idxq,1)
        for idxs,orbs in enumerate(sc):
            op_s = sqC(sc,idxs,0)
            op = op_q.dot(op_s)
            for idxp,orbp in enumerate(sl):
                for idxr,orbr in enumerate(sl):
                    tmp[idxp,idxr,0] += h2e[eint(orbp,orbq,orbr,orbs)]*op
    tmp = tmp.reshape((nl*nl,1,2**nc,2**nc))
    return tmp

#---------------------
# row-14 => T1[1,i-1]
#---------------------
# w[14,2]=v[pL,qR,rC,sC]*ta_rCsC
def l14r2(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nr,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxr,orbr in enumerate(sc):
        op_r = sqC(sc,idxr,0)
        for idxs,orbs in enumerate(sc):
            op_s = sqC(sc,idxs,0)
            if orbr < orbs:
                op = op_r.dot(op_s.dot(sgn))
                for idxp,orbp in enumerate(sl):
                    for idxq,orbq in enumerate(sr):
                        tmp[idxp,idxq] += h2e[eint(orbp,orbq,orbr,orbs)]*op
    return tmp

# w[14,3]=v[pL,qC,rC,sR]*(ta_qC^+rC)
def l14r3(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nr,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxq,orbq in enumerate(sc):
        op_q = sqC(sc,idxq,1)
        for idxr,orbr in enumerate(sc):
            op_r = sqC(sc,idxr,0)
            op = op_q.dot(op_r.dot(sgn))
            for idxp,orbp in enumerate(sl):
                for idxs,orbs in enumerate(sr):
                    tmp[idxp,idxs] += h2e[eint(orbp,orbq,orbr,orbs)]*op
    return tmp

# w[14,5]=v[pL,qC,rR,sR]*qC^+
def l14r5(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    if int(nr*(nr-1)/2)==0:
        print('error: mpo_consw.l14r5')
        exit()
    tmp = numpy.zeros((nl,int(nr*(nr-1)/2),2**nc,2**nc))
    for idxq,orbq in enumerate(sc):
        op_q = sqC(sc,idxq,1)
        for idxp,orbp in enumerate(sl):
            irs = 0 
            for idxr,orbr in enumerate(sr):
                for idxs,orbs in enumerate(sr):
                    if orbr < orbs:
                        tmp[idxp,irs] += h2e[eint(orbp,orbq,orbr,orbs)]*op_q
                        irs += 1
    return tmp

# w[14,9]=-a_rC
#---------------------------
# T1[p'] = M[p',pr]*Q[pr]
# M[p',pr]= d[p',p]*(-ar)
#---------------------------
def l14r9(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nl,nc,2**nc,2**nc))
    for idxr,orbr in enumerate(sc):
        op_r = sqC(sc,idxr,0)
        for idxp,orbp in enumerate(sl):
            tmp[idxp,idxp,idxr] = -op_r
    tmp = tmp.reshape((nl,nl*nc,2**nc,2**nc))
    return tmp

# w[14,12]=tI
def l14r12(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nl,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxq,orbq in enumerate(sl): 
        tmp[idxq,idxq] = sgn.copy()
    return tmp

# w[14,16]=T1
def l14r16(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,1,2**nc,2**nc))
    for idxq,orbq in enumerate(sc):
        op_q = sqC(sc,idxq,1)
        for idxr,orbr in enumerate(sc):
            op_r = sqC(sc,idxr,0)
            for idxs,orbs in enumerate(sc):
                op_s = sqC(sc,idxs,0)
                # r < s
                if orbr < orbs:
                    op = op_q.dot(op_r.dot(op_s))
                    for idxp,orbp in enumerate(sl):
                        tmp[idxp,0] += h2e[eint(orbp,orbq,orbr,orbs)]*op
    return tmp

#---------------------
# row-15 => T3[1,i-1]
#---------------------
# w[15,2]=-v[pC,qR,rL,sC]*ta_pC^+sC
def l15r2(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nr,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxp,orbp in enumerate(sc):
        op_p = sqC(sc,idxp,1)
        for idxs,orbs in enumerate(sc):
            op_s = sqC(sc,idxs,0)
            op = op_p.dot(op_s.dot(sgn))
            for idxr,orbr in enumerate(sl):
                for idxq,orbq in enumerate(sr):
                    tmp[idxr,idxq] += -h2e[eint(orbp,orbq,orbr,orbs)]*op
    return tmp

# w[15,3]=v[pC,qC,rL,sR]*(ta_pC^+qC^+)
def l15r3(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nr,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxq,orbq in enumerate(sc):
        op_q = sqC(sc,idxq,1)
        for idxp,orbp in enumerate(sc):
            op_p = sqC(sc,idxp,1)
            if orbp < orbq:
                op = op_p.dot(op_q.dot(sgn))
                for idxr,orbr in enumerate(sl):
                    for idxs,orbs in enumerate(sr):
                        tmp[idxr,idxs] += h2e[eint(orbp,orbq,orbr,orbs)]*op
    return tmp

# w[15,4]=v[pR,qR,rL,sC]*sC
def l15r4(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    if int(nr*(nr-1)/2)==0:
        print('error: mpo_consw.l15r4')
        exit()
    tmp = numpy.zeros((nl,int(nr*(nr-1)/2),2**nc,2**nc))
    for idxs,orbs in enumerate(sc):
        op_s = sqC(sc,idxs,0)
        for idxr,orbr in enumerate(sl):
            ipq = 0 
            for idxp,orbp in enumerate(sr):
                for idxq,orbq in enumerate(sr):
                    if orbp < orbq:
                        tmp[idxr,ipq] += h2e[eint(orbp,orbq,orbr,orbs)]*op_s
                        ipq += 1
    return tmp

# w[15,10]=a^+_p 
#-------------------------------
# T3[r'] = M[r',pr] * Q[pr]
# M[r',pr] = d[r',r](p^+_rC)
#-------------------------------
def l15r10(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nc,nl,2**nc,2**nc))
    for idxp,orbp in enumerate(sc):
        op_p = sqC(sc,idxp,1)
        for idxr,orbr in enumerate(sl):
            tmp[idxr,idxp,idxr] = op_p
    tmp = tmp.reshape((nl,nc*nl,2**nc,2**nc))
    return tmp

# w[15,14]=tI
def l15r14(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,nl,2**nc,2**nc))
    sgn = sgnC(sc)
    for idxq,orbq in enumerate(sl): 
        tmp[idxq,idxq] = sgn.copy()
    return tmp

# w[15,16]=T3
def l15r16(h1e,h2e,sl,sc,sr):
    nl = len(sl)
    nc = len(sc)
    nr = len(sr)
    tmp = numpy.zeros((nl,1,2**nc,2**nc))
    for idxp,orbp in enumerate(sc):
        op_p = sqC(sc,idxp,1)
        for idxq,orbq in enumerate(sc):
            op_q = sqC(sc,idxq,1)
            if orbp < orbq:
                for idxs,orbs in enumerate(sc):
                    op_s = sqC(sc,idxs,0)
                    op = op_p.dot(op_q.dot(op_s))
                    for idxr,orbr in enumerate(sl):
                         tmp[idxr,0] += h2e[eint(orbp,orbq,orbr,orbs)]*op
    return tmp

#------------------
# row-16 => I[i,K]
#------------------
# w[16,16]=I
def l16r16(h1e,h2e,sl,sc,sr):
    tmp = idnC(sc)
    return tmp
