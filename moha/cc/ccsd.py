from moha.cc.cc_auxiliary import *
from moha.io.log import log,timer

import numpy as np

__all__ = []

def ccsd(Eri,Ca,Cb,Emo_a,Emo_b,na,nb,I_thld,E_thld):
    """The second-order Moller-Plesset perturbation (MP2) algorithm
    implementation.

    Parameters
    ----------
    Eri : np.4darray
        Electron repulsion integral with spatial atomic orbitals.

    Ca : np.2darray
        Alpha molecular coefficient matrix.

    Cb : np.2darray
        Beta molecular coefficient matrix.

    Emo_a : np.1darray
        Alpha molecular orbital energy.

    Emo_b : np.1darray
        Beta molecular orbital energy.

    na : int
        Number of alpha electrons.

    nb : int
        Number of beta electrons.

    Returns
    -------
    results : dict
        Moller-Plesset perturbation results.
    """
    # Build the occupied and virtual orbitals slice
    norb = Eri.shape[0]
    oa = slice(0,2*na,2)
    ob = slice(1,2*nb,2)
    va = slice(2*na,2*norb,2)
    vb = slice(2*nb+1,2*norb,2)

    # Extend molecular orbital energy from spatial to spin basis
    fs = np.zeros(2*norb)
    fs[0::2] = Emo_a
    fs[1::2] = Emo_b

    # Transfer Eri tensor from atomic basis to molecular basis
    Eri_a = np.einsum("ijkl,ia,jb,kc,ld->abcd", Eri,Ca,Ca,Ca,Ca)
    Eri_b = np.einsum("ijkl,ia,jb,kc,ld->abcd", Eri,Cb,Cb,Cb,Cb)
    Eri_ab = np.einsum("ijkl,ia,jb,kc,ld->abcd", Eri,Ca,Ca,Cb,Cb)

    # Transfer Eri from spatial orbital to spin orbital
    Eri_a = np.repeat(Eri_a,2,axis=0)
    Eri_a = np.repeat(Eri_a,2,axis=1)
    Eri_a = np.repeat(Eri_a,2,axis=2)
    Eri_a = np.repeat(Eri_a,2,axis=3)
    node = np.zeros((2,2,2,2))
    node[0,0,0,0] = 1
    net = np.tile(node,(norb,norb,norb,norb))
    Eri = Eri_a*net

    Eri_b = np.repeat(Eri_b,2,axis=0)
    Eri_b = np.repeat(Eri_b,2,axis=1)
    Eri_b = np.repeat(Eri_b,2,axis=2)
    Eri_b = np.repeat(Eri_b,2,axis=3)
    node = np.zeros((2,2,2,2))
    node[1,1,1,1] = 1
    net = np.tile(node,(norb,norb,norb,norb))
    Eri += Eri_b*net

    Eri_ab = np.repeat(Eri_ab,2,axis=0)
    Eri_ab = np.repeat(Eri_ab,2,axis=1)
    Eri_ab = np.repeat(Eri_ab,2,axis=2)
    Eri_ab = np.repeat(Eri_ab,2,axis=3)
    node = np.zeros((2,2,2,2))
    node[0,0,1,1] = 1
    net = np.tile(node,(norb,norb,norb,norb))
    Eri += Eri_ab*net

    node = np.zeros((2,2,2,2))
    node[1,1,0,0] = 1
    net = np.tile(node,(norb,norb,norb,norb))
    Eri += Eri_ab.transpose(2,3,0,1)*net

    # Transfer Eri from chemists notation to physicists notation
    # <pq|rs> = [pr|qs]
    vs = Eri.transpose(0,2,1,3)

    # Slice Eri
    # <ab|cd>
    vs_vvvv = np.concatenate((vs[va,:,:,:],vs[vb,:,:,:]),axis=0)
    vs_vvvv = np.concatenate((vs_vvvv[:,va,:,:],vs_vvvv[:,vb,:,:]),axis=1)
    vs_vvvv = np.concatenate((vs_vvvv[:,:,va,:],vs_vvvv[:,:,vb,:]),axis=2)
    vs_vvvv = np.concatenate((vs_vvvv[:,:,:,va],vs_vvvv[:,:,:,vb]),axis=3)
    # <ib|cd>
    vs_ovvv = np.concatenate((vs[oa,:,:,:],vs[ob,:,:,:]),axis=0)
    vs_ovvv = np.concatenate((vs_ovvv[:,va,:,:],vs_ovvv[:,vb,:,:]),axis=1)
    vs_ovvv = np.concatenate((vs_ovvv[:,:,va,:],vs_ovvv[:,:,vb,:]),axis=2)
    vs_ovvv = np.concatenate((vs_ovvv[:,:,:,va],vs_ovvv[:,:,:,vb]),axis=3)
    # <ij|cd>
    vs_oovv = np.concatenate((vs[oa,:,:,:],vs[ob,:,:,:]),axis=0)
    vs_oovv = np.concatenate((vs_oovv[:,oa,:,:],vs_oovv[:,ob,:,:]),axis=1)
    vs_oovv = np.concatenate((vs_oovv[:,:,va,:],vs_oovv[:,:,vb,:]),axis=2)
    vs_oovv = np.concatenate((vs_oovv[:,:,:,va],vs_oovv[:,:,:,vb]),axis=3)
    # <ij|kd>
    vs_ooov = np.concatenate((vs[oa,:,:,:],vs[ob,:,:,:]),axis=0)
    vs_ooov = np.concatenate((vs_ooov[:,oa,:,:],vs_ooov[:,ob,:,:]),axis=1)
    vs_ooov = np.concatenate((vs_ooov[:,:,oa,:],vs_ooov[:,:,ob,:]),axis=2)
    vs_ooov = np.concatenate((vs_ooov[:,:,:,va],vs_ooov[:,:,:,vb]),axis=3)
    # <ij|kl>
    vs_oooo = np.concatenate((vs[oa,:,:,:],vs[ob,:,:,:]),axis=0)
    vs_oooo = np.concatenate((vs_oooo[:,oa,:,:],vs_oooo[:,ob,:,:]),axis=1)
    vs_oooo = np.concatenate((vs_oooo[:,:,oa,:],vs_oooo[:,:,ob,:]),axis=2)
    vs_oooo = np.concatenate((vs_oooo[:,:,:,oa],vs_oooo[:,:,:,ob]),axis=3)
    
    # Build the orbital energy denominator
    fso = np.concatenate((fs[oa],fs[ob]))
    fsv = np.concatenate((fs[va],fs[vb]))
    Ds = fso.reshape(-1,1)-fsv
    Dd = fso.reshape(-1,1,1,1)+fso.reshape(-1,1,1)-fsv.reshape(-1,1)-fsv

    # Initial guess
    ts = np.zeros((na+nb, 2*norb-na-nb))
    td = vs_oovv/Dd
    Eold = energy(fso,fsv,vs_oovv,ts,td)

    ### Start CCSD iterations
    for Iter in range(1,I_thld+1):
        # Build intermediates
        Ivo = intermediate_Ivo(fso,fsv,vs_oovv,ts)
        Ivv = intermediate_Ivv(fso,fsv,vs_ovvv,vs_oovv,ts,td)
        IPoo = intermediate_IPoo(fso,vs_oovv,vs_ooov,ts,td)
        Ioo = intermediate_Ioo(IPoo,Ivo,ts)
        Ioooo = intermediate_Ioooo(vs_oooo,vs_oovv,vs_ooov,ts,td)
        Iovov = intermediate_Iovov(vs_ovvv,vs_oovv,vs_ooov,ts,td)
        IPvvvo = intermediate_IPvvvo(vs_ovvv,vs_oovv,ts)
        IPovoo = intermediate_IPovoo(vs_ovvv,vs_ooov,ts,td)
        coovv = intermediate_coovv(ts,td)
        
        # Update 
        ts = update_ts(fso,fsv,vs_ovvv,vs_oovv,vs_ooov,ts,td,Ivv,IPoo,Ivo)
        td = update_td(fso,fsv,vs_oovv,vs_vvvv,ts,td,Ivv,Ioo,Ioooo,Iovov,IPvvvo,IPovoo)

        # Update energy
        E = energy(fso,fsv,vs_oovv,ts,td)

        # Update convergence
        dE = abs(E-Eold)

        # Store Updated energy
        Eold = E

        log('{0:<20d} {1:<20.6f} {2:<20.6f}'.format(Iter, E, dE))
        if dE < E_thld:
            break
        if Iter == I_thld:
            raise Exception("Maximum number of SCF cycles exceeded.")

    return E
