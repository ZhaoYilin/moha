import numpy as np

################################################################################
# Three F intermediates: Fvv,Foo,Fov                                    
# Three W intermediates: Woooo,Wvvvv,Wovvo
################################################################################
def intermediate_Fvv(fso,fsv,vs_ovvv,vs_oovv,ts,td):
    """
    Parameters
    ----------
    fs : np.1darray
        Fock matrix diagonal with spin orbitals.

    vs : np.4darray
        Electron repulsion integral with spin orbitals.
    
    ts : np.2darray
        Single cluster amplitudes.
    """
    fvv = fsv.reshape(-1,1)-fsv
    fvv[np.diag_indices_from(fvv)] = 0
    fov = fsv.reshape(-1,1)-fso
    ttau = effective_tilde_tau(ts,td)

    Fvv = fvv
    Fvv -= 0.5*np.einsum('me,ma->ae', fov,ts)
    Fvv += np.einsum('mf,mafe->ae', ts,vs_ovvv)
    Fvv -= 0.5*np.einsum('mnaf,mnef->ae', ttau,vs_oovv)

    return Fvv

def intermediate_Foo(fso,fsv,vs_oovv,vs_ooov,ts,td):
    """
    Parameters
    ----------
    f : np.1darray
        fock spin orbital

    v : np.4darray
        fock spin orbital
    
    ts : np.2darray
        fock spin orbital

    td : np.4darray
        fock spin orbital
    """
    foo = fso.reshape(-1,1)-fso
    foo[np.diag_indices_from(foo)] = 0
    fov = fsv.reshape(-1,1)-fso
    ttau = effective_tilde_tau(ts,td)

    Foo = foo
    Foo += 0.5*np.einsum('ie,me->mi', ts,fov)
    Foo += np.einsum('ne,mnie->mi', ts,vs_ooov)
    Foo += 0.5*np.einsum('inef,mnef->mi', ttau,vs_oovv)

    return Foo


def intermediate_Fov(fso,fsv,vs_oovv,ts):
    """
    Parameters
    ----------
    f : np.1darray
        fock spin orbital

    v : np.4darray
        fock spin orbital
    
    ts : np.2darray
        fock spin orbital
    """
    fov = fsv.reshape(-1,1)-fso

    Fov = fov
    Fov += np.einsum('nf,mnef->me', ts,vs_oovv)

    return Fov

def intermediate_Woooo(vs_oovv,vs_ooov,vs_oooo,ts,td):
    tmp = np.einsum('je,mnie->mnij', ts,vs_ooov)
    tau = effective_tau(ts,td)

    Woooo = vs_oooo
    Woooo += tmp
    Woooo -= tmp.transpose(0,1,3,2)
    Woooo += 0.25*np.einsum('ijef,mnef->mnij', tau,vs_oovv)

    return Woooo

def intermediate_Wvvvv(vs_vvvv,vs_oovv,vs_ooov,ts,td):
    tmp = np.einsum('je,mnie->mnij', ts,vs_ooov)
    tau = effective_tau(ts,td)

    Wvvvv = vs_vvvv
    Wvvvv -= tmp
    Wvvvv += tmp.transpose(1,0,2,3)
    Wvvvv += 0.25*np.einsum('mnab,mnef->abef', tau,vs_oovv)

    return Wvvvv

def intermediate_Wovvo(vs_ovvv,vs_oovv,vs_ooov,ts,td):
    vs_ovvo = vs_oovv.transpose(0,3,2,1)
    vs_oovo = vs_ooov.transpose(0,1,3,2)

    Wovvo = vs_ovvo
    Wovvo += np.einsum('jf,mbef->mbej', ts,vs_ovvv)
    Wovvo -= np.einsum('nb,mnej->mbej', ts,vs_oovo)
    Wovvo -= 0.5*np.einsum('jnfb,mnef->mbej', td,vs_oovv)
    Wovvo -= np.einsum('jf,nb,mnej->mbej', ts,ts,vs_oovv)

    return Wovvo


################################################################################
# Effective two-particle excitation operators
################################################################################
def effective_tilde_tau(ts,td):
    """
    Parameters
    ----------
    ts : np.1darray
        fock spin orbital

    td : np.4darray
    """
    ttau = td.copy()
    ttau += 0.5*np.einsum('ia,jb->ijab',ts,ts)
    ttau -= 0.5*np.einsum('ib,ja->ijab',ts,ts)
    return ttau

def effective_tau(ts,td):
    """
    """
    tau = td.copy()
    tau += np.einsum('ia,jb->ijab',ts,ts)
    tau -= np.einsum('ib,ja->ijab',ts,ts)
    return tau

################################################################################
# MakeT1 and makeT2, as they imply, construct the actual amplitudes
################################################################################
def amplitude_ts(fso,fsv,vs_ovvv,vs_oovv,vs_ooov,ts,td):
    """amplitudes
    Ref Stanton eq(1)
    """
    Fvv = intermediate_Fvv(fso,fsv,vs_ovvv,vs_oovv,ts)
    Foo = intermediate_Foo(fso,fsv,vs_oovv,vs_ooov,ts,td)
    Fov = intermediate_Fov(fso,fsv,vs_oovv,ts)

    Dsts = fso.reshape(-1,1) - fsv
    Dsts += np.einsum('ie,ae->ia',ts,Fvv)
    Dsts -= np.einsum('ma,mi->ia',ts,Foo)
    Dsts += np.einsum('imae,me->ia',td,Fov)
    Dsts -= np.einsum('nf,naif->ia',ts,vs_oovv.transpose(0,2,1,3)) #ovov
    Dsts -= 0.5*np.einsum('imef,maef->ia',td,vs_ovvv)
    Dsts -= 0.5*np.einsum('mnae,nmei->ia',td,vs_ooov.transpose(0,1,3,2)) #vs_oovo

    Ds = fso.reshape(-1,1) - fsv
    ts = Dsts/Ds

    return ts

def amplitude_td(fso,fsv,vs_vvvv,vs_ovvv,vs_oovv,vs_ooov,vs_oooo,ts,td):
    """amplitudes
    Ref Stanton eq(1)
    """
    Fvv = intermediate_Fvv(fso,fsv,vs_ovvv,vs_oovv,ts)
    Foo = intermediate_Foo(fso,fsv,vs_oovv,vs_ooov,ts,td)
    Fov = intermediate_Fov(fso,fsv,vs_oovv,ts)
    Woooo = intermediate_Woooo(vs_oovv,vs_ooov,vs_oooo,ts,td)
    Wvvvv = intermediate_Wvvvv(vs_vvvv,vs_oovv,vs_ooov,ts,td)
    Wovvo = intermediate_Wovvo(vs_ovvv,vs_oovv,vs_ooov,ts,td)


    Ddtd = vs_oovv

    # P_(ab)
    tmp1 = np.einsum('ijae,be->ijab',td,Fvv)
    tmp1 -= 0.5*np.einsum('ijae,mb,me->ijab',td,ts,Fov)
    Ddtd += tmp1
    Ddtd -= tmp1.transpose(0,1,3,2)

    # P_(ij)
    tmp2 = np.einsum('imab,mj->ijab',td,Foo)
    tmp2 += 0.5*np.einsum('imab,je,me->ijab',td,ts,Fov)
    Ddtd += tmp2
    Ddtd -= tmp2.transpose(1,0,2,3)

    Ddtd += 0.5*np.einsum('mnab,mnij->ijab',tau,Woooo)
    Ddtd += 0.5*np.einsum('ijef,abef->ijab',tau,Wvvvv)
    
    # P_(ij) * P_(ab)
    tmp3 = np.einsum('imae,mbej->ijab',td,Wovvo)
    tmp3 -= np.einsum('ie,ma,mbej->ijab',ts,ts,vs_ovvo)
    Ddtd += tmp3 #ijab
    Ddtd -= tmp3.transpose(0,1,3,2) #ijba
    Ddtd -= tmp3.transpose(1,0,2,3) #jiab
    Ddtd += tmp3.transpose(1,0,3,2) #jiba

    # P_(ij)
    tmp4 = np.einsum('ie,abej->ijab',ts,vs_vvvo)
    Ddtd += tmp4
    Ddtd -= tmp4.transpose(1,0,2,3)

    # P_(ab)
    tmp5 = np.einsum('ma,mbij->ijab',ts,vs_ovoo)
    Ddtd += tmp5
    Ddtd -= tmp5.transpose(0,1,3,2)

    Dd = fso.reshape(-1,1,1,1)+fso.reshape(-1,1,1)-fsv.reshape(-1,1)-fsv
    td = Ddtd/Dd

    return td

################################################################################
# Compute the CCSD energy
################################################################################
def energy(fso,fsv,vs_oovv,ts,td):
    ### Compute CCSD correlation energy
    double_bar = vs_oovv - vs_oovv.transpose(0,1,3,2)
    fov = fso.reshape(-1,1)-fsv
    E = np.einsum('ia,ia->', fov, ts)
    E += 0.5 * np.einsum('ijab,ia,jb->', double_bar, ts, ts)
    print('E',E)
    E += 0.25 * np.einsum('ijab,ijab->', double_bar, td)
    return E
