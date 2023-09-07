import itertools
import sys

from moha.symmetry.point_group.locate_symmetry_operator import *
from moha.symmetry.point_group.symmetry_operator import *
import numpy as np

__all__ = ["PointGroup"]


class PointGroup(set):
    """ A special case of group based on a set of symmetry operator {E,i,sigma,Cn,Sn...},
        together with an dot operation, satisfying the following properties:
    
        1. Closedness
        2. Associativity
        3. Identity element
        4. Inverses

    Attributes
    ----------
    label : str
        Standard Schönflies symbol. 

    Properties
    ----------
    E : List
        List of identity symmetry operators.

    i : List
        List of inversion symmetry operators.

    sigma : List
        List of reflection symmetry operators.

    Cn : List
        List of rotation symmetry operators.

    Sn : List
        List of improper rotation symmetry operators.

    Cn_principals : List
        List of rotation operators with the highest order n.

    C2_perpendicular : List
        List of proper rotation operators with order 2 which is perpendicular to the 
        principal axis.

    sigma_h : List
        List of horizontal refelction symmetry operators which is perpendicular to the 
        principal axis.

    sigma_v : List
        List of vertical reflection symmetry operators which includes the principal axis
        and passes through the bonds.

    sigma_d : List
        List of dihedral reflection symmetry operators which is vertical reflection 
        symmetry operators which includes the principal axis while bisecting the angle 
        between two C2 axes that are perpendicular to it. 

    cardinality : int
        Number of symmetry operators in the point group.

    irreducible_cardinality : int
        Number of irreducible representations.

    is_group : bool
        Check if the instance is a valid group.

    is_abelian_group : bool
        Check if the instance is a valid abelian group.

    subgroups : List
        List of subgroups for the group.

    D2h_subgroup : PointGroup
        The highest group within D2h and subgroups of D2h.

    Methods
    -------
    build_label(self)
        Generate the label of the point group by analysis the symmetry elements.

    dot(self, operator_a, operator_b)
        The binary operation that can be applied to two elements of the set.

    ClassMethods
    ------------
    build(cls, mol, error=1e-3)
        Build the poin group for a given molecule.
    """
    def __new__(cls, symmetry_operators):
        """ Generate new molecule object.

        Parameters
        ----------
        symmetry_operators : List
            List of symmetry operators.
        """
        assert all(
            isinstance(operator, SymmetryOperator)
            for operator in symmetry_operators)
        obj = set().__new__(cls, symmetry_operators)

        return obj

    def __init__(self, symmetry_operators, label=None):
        """ Initialize the instance.

        Parameters
        ----------
        symmetry_operators : List
            List of symmetry operators.
        """
        super().__init__(symmetry_operators)
        if label == None:
            self.label = self.build_label()
        elif isinstance(label, str):
            self.label = label
        else:
            raise TypeError("Label should be None or str.")

    @property
    def E(self):
        """ The proper rotation operators with the highest order n.

        Return
        ------
        Cn_list : List
            List of principal proper rotaiton operators.
        """
        operator_list = []
        for operator in self:
            if isinstance(operator, Identity):
                operator_list.append(operator)
        return operator_list

    @property
    def i(self):
        """ The proper rotation operators with the highest order n.

        Return
        ------
        Cn_list : List
            List of principal proper rotaiton operators.
        """
        operator_list = []
        for operator in self:
            if isinstance(operator, Inversion):
                operator_list.append(operator)
        return operator_list

    @property
    def sigma(self):
        """ The proper rotation operators with the highest order n.

        Return
        ------
        Cn_list : List
            List of principal proper rotaiton operators.
        """
        operator_list = []
        for operator in self:
            if isinstance(operator, Reflection):
                operator_list.append(operator)
        return operator_list

    @property
    def Cn(self):
        """ The proper rotation operators with the highest order n.

        Return
        ------
        Cn_list : List
            List of principal proper rotaiton operators.
        """
        operator_list = []
        for operator in self:
            if isinstance(operator, Rotation):
                operator_list.append(operator)
        return operator_list

    @property
    def Sn(self):
        """ The proper rotation operators with the highest order n.

        Return
        ------
        Cn_list : List
            List of principal proper rotaiton operators.
        """
        operator_list = []
        for operator in self:
            if isinstance(operator, ImproperRotation):
                operator_list.append(operator)
        return operator_list

    @property
    def Cn_principals(self):
        """ The proper rotation operators with the highest order n.

        Return
        ------
        Cn_list : List
            List of principal proper rotaiton operators.
        """
        if len(self.Cn) > 0:
            order_list = [Cn.order for Cn in self.Cn]
            Cn_indice = np.where(np.array(order_list) == max(order_list))
            Cn_list = np.array(self.Cn)[Cn_indice]
            return Cn_list
        else:
            return []

    @property
    def C2_perpendicular(self):
        """ The proper rotation operators with order 2 which is perpendicular to the 
            principal axis.

        Return
        ------
        C2_list : List
            List of order 2 proper rotaiton operators.
        """
        C2_perpendicular_list = []
        if len(self.Cn) > 0:
            order_list = [Cn.order for Cn in self.Cn]
            C2_indice = np.where(np.array(order_list) == 2)
            C2_list = np.array(self.Cn)[C2_indice]

            principal_axis = self.Cn_principals[0].symmetry_element
            for C2 in C2_list:
                norm_vector = C2.symmetry_element
                # Is orthogonal
                if np.isclose(np.dot(principal_axis, norm_vector), 0):
                    C2_perpendicular_list.append(C2)
        return C2_perpendicular_list

    @property
    def sigma_h(self):
        """ Horizontal refelction plane which is perpendicular to the principal axis.
        """
        sigma_h_list = []
        if len(self.sigma) > 0:
            principal_axis = self.Cn_principals[0].symmetry_element
            for sigma in self.sigma:
                norm_vector = sigma.symmetry_element
                # Is parallel
                if np.allclose(np.cross(principal_axis, norm_vector),
                               [0, 0, 0]):
                    sigma_h_list.append(sigma)
        return sigma_h_list

    @property
    def sigma_v(self):
        """ Vertical reflection symmetry operators which includes the principal axis and passes 
        through the bonds.

        Return
        ------
        C2_list : List
            List of order 2 proper rotaiton operators.
        """
        sigma_v_list = []
        if len(self.sigma) > 0:
            principal_axis = self.Cn_principals[0].symmetry_element
            for sigma in self.sigma:
                norm_vector = sigma.symmetry_element
                # Is orthogonal
                if np.isclose(np.dot(principal_axis, norm_vector),
                              0,
                              rtol=1e-3):
                    sigma_v_list.append(sigma)
        return sigma_v_list

    @property
    def sigma_d(self):
        """ Dihedral reflection plane is a special case of vertical reflection plane 
        which includes the principal axis while bisecting the angle between two C2 axes
        that are perpendicular to it. 
        """
        sigma_d_list = []
        if len(self.sigma) > 0:
            C2_perpendicular_list = self.C2_perpendicular
            sigma_v_list = self.sigma_v
            for sigma in sigma_v_list:
                for C2_a, C2_b in itertools.combinations(
                        C2_perpendicular_list, 2):
                    axis_a = C2_a.symmetry_element
                    axis_b = C2_b.symmetry_element
                    norm_vector = sigma.symmetry_element
                    angle_na = angle(axis_a, norm_vector)
                    angle_nb = angle(axis_b, norm_vector)
                    # Is orthogonal
                    if np.isclose(angle_na, angle_nb, rtol=1e-3):
                        sigma_d_list.append(sigma)
        return sigma_d_list

    @property
    def order(self):
        """ Number of symmetry operators present in the point group.

        Returns
        -------
        n : int
            Number of the symmetry operators.
        """
        return len(self)

    @property
    def rank(self):
        """ Number of symmetry operators present in the generator subgroup.

        Returns
        -------
        n : int
            Number of the symmetry operators.
        """
        return len(self)

    @property
    def irreducible_cardinality(self):
        """ Number of irreducible representations.

        Returns
        -------
        n : int
            Number of the irreducible symmetry operators.
        """
        return len(self)

    @property
    def is_group(self):
        """ Check if the instance is a valid group satisfying the following properties:

            1. Closedness
            2. Associativity
            3. Identity element
            4. Inverses
        """
        # Identity element
        if not Identity() in self:
            return False
        # Closedness
        for a, b in itertools.product(self, repeat=2):
            if not self.dot(a, b) in self:
                return False
        # Associativity
        for a, b, c in itertools.product(self, repeat=3):
            if not self.dot(self.dot(a, b), c) == self.dot(a, self.dot(b, c)):
                return False
        # Inverses
        for a in self:
            if not Identity() in [self.dot(a, b) for b in self]:
                return False
        else:
            return True

    @property
    def is_abelian_group(self):
        """ Check if the instance is a valid abelian group satisfying the following 
            properties:

            1. Closedness
            2. Associativity
            3. Identity element
            4. Inverses
            5. Commutative
        """
        if self.is_group:
            for a, b in itertools.product(self, repeat=2):
                if self.dot(a, b) == self.dot(b, a):
                    return True
                else:
                    return False
        else:
            return False

    @property
    def subgroups(self):
        """ List of subgroups for the group.
        """
        subgroups = []
        subsets = []
        for n in range(len(self) + 1)[1:]:
            for subset in itertools.combinations(self, n):
                subsets.append(set(subset))
        for subset in subsets:
            subgroup_candidate = self.__class__(subset)
            if subgroup_candidate.is_group:
                subgroups.append(subgroup_candidate)

        return subgroups

    @property
    def D2h_subgroups(self):
        """ The highest group within D2h and subgroups of D2h.
        """
        two_fold_operators = []
        D2h_subgroups = []
        for operator in self:
            if self.dot(operator, operator) == Identity():
                two_fold_operators.append(operator)
        group = self.__class__(two_fold_operators)
        for subgroup in group.subgroups:
            if subgroup.is_abelian_group:
                D2h_subgroups.append(subgroup)

        return D2h_subgroups

    @classmethod
    def build(cls, mol, error=1e-3):
        """ Build the poin group for a given molecule.

        Parameters
        ----------
        mol : Molecule
            The molecule instance.

        Returns
        -------
        group : PointGroup
            The pointgroup for the given molecule.

        mol : Molecule
            The molecule whose center of mass move to origin with standary orientation.
        """
        # Move the center of mass to origin.
        mol = mol.translation(mol.center_of_mass)
        # Set the standary orientation.
        mol = mol.standary_orientation(error)

        symmetry_operators = []
        if not mol.pg:
            symmetry_operators.append(Identity())
            group = cls(symmetry_operators, "C_{1}")
            return group, mol

        elif len(mol) == 1:
            symmetry_operators.append(Identity())
            group = cls(symmetry_operators, "SO(3)")
            return group, mol

        # Linear
        # Fake a highest subgoup within D2h.
        error_order = -int(math.log10(error))
        v, w = mol.principal_moments_of_inertia
        v_round = np.around(v, error_order)
        if v[0] == 0 and v[1] == v[2]:
            if Inversion()(mol) == mol:
                E_list = [Identity()]
                i_list = [Inversion()]
                sigma_list = [
                    Reflection([1.0, 0.0, 0.0]),
                    Reflection([0.0, 1.0, 0.0]),
                    Reflection([0.0, 0.0, 1.0])
                ]
                Cn_list = [
                    Rotation([1.0, 0.0, 0.0], 2),
                    Rotation([0.0, 1.0, 0.0], 2),
                    Rotation([0.0, 0.0, 1.0], 2)
                ]
                group = cls(E_list + i_list + sigma_list + Cn_list, "D_{ooh}")
                return group, mol
            else:
                E_list = [Identity()]
                i_list = [Inversion()]
                sigma_list = [Reflection([0.0, 0.0, 1.0])]
                Cn_list = [Rotation([0.0, 0.0, 1.0], 2)]
                group = cls(E_list + i_list + sigma_list + Cn_list, "C_{oov}")
                return group, mol

        else:
            # Locating symmetry operators
            # identity_operators
            symmetry_operators += [Identity()]
            # inversion_operators
            if Inversion()(mol) == mol:
                symmetry_operators += [Inversion()]
            #rotation_operators
            rotation_operators = locate_rotation(mol, error)
            for i, Cn in enumerate(rotation_operators):
                rotation_operators[i] = Cn.expansion
            rotation_operators = list(itertools.chain(*rotation_operators))
            symmetry_operators += rotation_operators
            #improper_rotation_operators
            improper_rotation_operators = locate_improper_rotation(
                mol, rotation_operators)
            symmetry_operators += improper_rotation_operators
            # reflection_operators
            reflection_operators = locate_reflection(mol, error)
            symmetry_operators += reflection_operators

            group = cls(symmetry_operators)
            return group, mol

    def build_label(self):
        """ Generate the label of the point group by analysis the symmetry elements.

        Return
        ------
        label : str
            Standard Schönflies symbol. 
        """
        # Empty
        if len(self) == 0:
            return "Empty"

        # Check if it is a valid group.
        #if not self.is_group:
        #    raise AssertionError("The given set is not a valid group.")

        # Low symmetry
        if len(self.Cn) == 0:
            if len(self.sigma) == 1:
                return "C_{s}"
            elif len(self.i) == 1:
                return "C_{i}"
            else:
                return "C_{1}"

        # High symmetry
        if len(self.Cn_principals) > 2 and self.Cn_principals[0].order > 2:
            if not len(self.i) == 1:
                return "T_{d}"
            elif any([Cn.order == 5 for Cn in self.Cn_prinpicals]):
                return "T_{h}"
            else:
                return "O_{h}"

        n = self.Cn_principals[0].order

        # Dihedral
        if len(self.C2_perpendicular) > 1:
            if len(self.sigma_h) > 0:
                return "D_{" + str(n) + "h}"
            elif len(self.sigma_v) > 0:
                return "D_{" + str(n) + "d}"
            else:
                return "D_{" + str(n) + "}"

        # Cyclic
        else:
            if len(self.sigma_h) > 0:
                return "C_{" + str(n) + "h}"
            elif len(self.sigma_v) > 0:
                return "C_{" + str(n) + "v}"
            elif len(self.Sn) > 0:
                return "S_{2" + str(n) + "}"
            else:
                return "C_{" + str(n) + "}"

    def dot(self, operator_a, operator_b):
        """ The binary operation that can be applied to two elements of the set.

        Parameters
        ----------
        operator_a : SymmetryOperator
            The element in the set.

        operator_b : SymmetryOperator
            The element in the set.

        Retrun
        ------
        operator : SymmetryOperator
            The element in the set.
        """
        matrix = np.dot(operator_a.matrix, operator_b.matrix)
        for operator in self:
            if np.allclose(operator.matrix, matrix, rtol=1e-3):
                return operator
