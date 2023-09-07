import itertools
import math

from moha.symmetry.point_group.auxiliary import *
from moha.symmetry.point_group.symmetry_operator import *
import numpy as np


def locate_Cn(sea, error=1e-3):
    """ Locate the proper rotations for each set of SEA based on the analysis of their
    number of atoms k, and the principal moments of inertia of the set.

    Parameters
    ----------
    sea : Molecule
        Subsets of symmetrically equivalent atoms.

    threshold : float
        Threshold of nearly equal.
    
    Return
    ------
    Cn_list : List
        List of the proper rotation symmetry elements.
    """
    k = len(sea)
    Cn_list = []
    error_order = -int(math.log10(error))
    v, w = sea.principal_moments_of_inertia
    v_round = np.around(v, error_order)
    if k == 1:
        # Single atom arrangement: Trivial C2 cross this atom
        return Cn_list

    elif k == 2:
        # Linear arrangement: Any order C2 and C2_p
        C_ = Rotation(w[:, 0], 2)
        Cn_list.append(C_)
        C_ = Rotation(w[:, 1], 2)
        Cn_list.append(C_)
        C_ = Rotation(w[:, 2], 2)
        print('wwfwf', w[:, 2])
        Cn_list.append(C_)

    elif k > 2:
        if np.isclose(v[0] + v[1], v[2], rtol=error):
            # Polygonal arrangement
            if v_round[0] == v_round[1]:
                # Regular polygonal: Ck and its divisor
                print("Regular polygonal: Ck and its divisor")
                C_ = Rotation(w[:, 2], k)
                Cn_list.append(C_)
                for n in divisors(k):
                    C_ = Rotation(w[:, 2], n)
                    Cn_list.append(C_)
            elif v_round[0] != v_round[1] != v_round[2]:
                # Ia != Ib != Ic
                # Irregular: Divisor of Ck
                print("Irregular: Divisor of Ck")
                for n in divisors(k):
                    C_ = Rotation(w[:, 2], n)
                    Cn_list.append(C_)
            else:
                return False
        else:
            # Polyhedral Arrangement
            print("Polyhedral Arrangement")
            if v_round[0] != v_round[1] != v_round[2]:
                # C2
                C_ = Rotation(w[:, 1], 2)
                Cn_list.append(C_)
            elif v_round[0] == v_round[1] < v_round[2]:
                # Ia = Ib < Ic
                # C_{k/2} and its divisor
                C_ = Rotation(w[:, 2], int(k / 2))
                Cn_list.append(C_)
                for n in divisors(int(k / 2)):
                    C_ = Rotation(w[:, 2], n)
                    Cn_list.append(C_)
            elif v_round[0] < v_round[1] == v_round[2]:
                # Ia < Ib = Ic
                # C_{k/2} and its divisor
                C_ = Rotation(w[:, 2], int(k / 2))
                Cn_list.append(C_)
                for n in divisors(int(k / 2)):
                    C_ = Rotation(w[:, 0], n)
                    Cn_list.append(C_)
            else:
                return False
    else:
        return False
    return Cn_list


def locate_C2(mol, error=1e-3):
    """ Locate the proper rotations for each set of SEA based on the analysis of their
    number of atoms k, and the principal moments of inertia of the set.

    Parameters
    ----------
    sea : Molecule
        Subsets of symmetrically equivalent atoms.

    threshold : float
        Threshold of nearly equal.
    
    Return
    ------
    Cn_list : List
        List of the proper rotation symmetry elements.
    """
    C2_list = []
    SEA = mol.symmetrically_equivalent_atoms()

    same_point_list = []
    for sea in SEA[1:]:
        same_point = np.allclose(SEA[0].center_of_mass,
                                 sea.center_of_mass,
                                 atol=error)
        same_point_list.append(same_point)
    if all(same_point_list):
        for sea in SEA:
            for i, j in itertools.product(range(len(sea)), repeat=2):
                if i < j:
                    atom_i, atom_j = sea[i], sea[j]
                    r_ij = (np.array(atom_i.coordinate) +
                            np.array(atom_j.coordinate)) / 2.0
                    r_CM = sea.center_of_mass
                    r_C2 = r_ij - r_CM
                    C2_list.append(Rotation(r_C2, 2))
            for atom_i in sea:
                r_i = np.array(atom_i.coordinate)
                r_CM = sea.center_of_mass
                r_C2 = r_i - r_CM
                C2_list.append(Rotation(r_C2, 2))

    return C2_list


def locate_principal_rotation(mol, error=1e-3):
    """ Locate principal rotation.
    """
    SEA = mol.symmetrically_equivalent_atoms()
    Cn_candidates = []
    Cn_verified = []
    for sea in SEA:
        Cn_candidates += locate_Cn(sea, error=1e-3)
    for Cn in Cn_candidates:
        Cn_verified.append(Cn)
    order_list = [Cn.order for Cn in Cn_verified]
    Cn_principals_indice = np.where(np.array(order_list) == max(order_list))
    Cn_principals = np.array(Cn_verified)[Cn_principals_indice]
    return list(Cn_principals)


def locate_rotation(mol, error):
    Cp_list = locate_principal_rotation(mol, error=1e-3)
    C2_list = locate_C2(mol, error=1e-3)
    rotation_list = Cp_list + C2_list
    return rotation_list


def locate_improper_rotation(mol, Cn_list):
    """
    """
    improper_rotation_list = []
    Sn_candidates = []
    for Cn in Cn_list:
        Sn = ImproperRotation(Cn.symmetry_element, Cn.order)
        S2n = ImproperRotation(Cn.symmetry_element, Cn.order * 2)
        Sn_candidates += [Sn, S2n]
    for Sn in Sn_candidates:
        if Sn(mol) == mol:
            improper_rotation_list.append(Sn)

    return improper_rotation_list


def locate_reflection(mol, error):
    """ locate proper rotations in sets of SEA based on the analysis of their
     number, k, and principal moments of inertia,
    """
    reflection_list = []
    sigma_candidates = []
    SEA = mol.symmetrically_equivalent_atoms()
    for sea in SEA:
        for i, j in itertools.product(range(len(sea)), repeat=2):
            if i < j:
                atom_i, atom_j = sea[i], sea[j]
                r_ij = (np.array(atom_i.coordinate) -
                        np.array(atom_j.coordinate))
                n = r_ij / np.linalg.norm(r_ij)
                sigma_candidates.append(Reflection(n))
    for sigma in sigma_candidates:
        if sigma(mol) == mol:
            reflection_list.append(sigma)

    return reflection_list
