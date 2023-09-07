import numpy as np

__all__ = []

def mp2(Eri,Ca,Cb,Emo_a,Emo_b,na,nb):
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
    # Build double bar 
    # <pq||rs> = <pq|rs> - <pq|sr>
    double_bar = Eri.transpose(0,2,1,3) - Eri.transpose(0,2,3,1)

    # Slice double bar <pq||rs>
    # <pq||rs> -> <ij||ab>
    db_oovv = np.concatenate((double_bar[oa,:,:,:],double_bar[ob,:,:,:]),axis=0)
    db_oovv = np.concatenate((db_oovv[:,oa,:,:],db_oovv[:,ob,:,:]),axis=1)
    db_oovv = np.concatenate((db_oovv[:,:,va,:],db_oovv[:,:,vb,:]),axis=2)
    db_oovv = np.concatenate((db_oovv[:,:,:,va],db_oovv[:,:,:,vb]),axis=3)
    # <pq||rs> -> <ij||ij>
    db_oooo = np.concatenate((double_bar[oa,:,:,:],double_bar[ob,:,:,:]),axis=0)
    db_oooo = np.concatenate((db_oooo[:,oa,:,:],db_oooo[:,ob,:,:]),axis=1)
    db_oooo = np.concatenate((db_oooo[:,:,oa,:],db_oooo[:,:,ob,:]),axis=2)
    db_oooo = np.concatenate((db_oooo[:,:,:,oa],db_oooo[:,:,:,ob]),axis=3)
    
    # Build the orbital energy denominator
    fso = np.concatenate((fs[oa],fs[ob]))
    fsv = np.concatenate((fs[va],fs[vb]))
    denom = fso.reshape(-1,1,1,1)+fso.reshape(-1,1,1)-fsv.reshape(-1,1)-fsv

    # Compute the energy
    Emp0 = np.einsum('i->',fso)
    Emp1 = -0.5*np.einsum('ijij->',db_oooo)
    Emp2 = 0.25*np.einsum('ijab,ijab->',db_oovv,db_oovv/denom)

    results = {
    "success": True,
    "Emp0": Emp0,
    "Emp1": Emp1,
    "Emp2": Emp2
    }
    return results
