from moha.symmetry.point_group.point_group import PointGroup
from moha.symmetry.point_group.symmetry_operator import *

class D2h(PointGroup):

    def __init__(self):
        identity_list = [
            Identity()
        ]

        inversion_list = [
            Inversion()
        ]

        reflection_list = [
            Reflection([1.0,0.0,0.0]),
            Reflection([0.0,1.0,0.0]),
            Reflection([0.0,0.0,1.0])
        ]

        rotation_list = [
            Rotation([1.0,0.0,0.0],2),
            Rotation([0.0,1.0,0.0],2),
            Rotation([0.0,0.0,1.0],2)
        ]


        improper_rotation_list = [
        ]

        label = 'D_{2h}'

        super().__init__(identity_list, inversion_list, reflection_list, rotation_list, improper_rotation_list, label)


class C2h(PointGroup):

    def __init__(self):
        identity_list = [
            Identity()
        ]

        inversion_list = [
            Inversion()
        ]

        reflection_list = [
            Reflection([0.0,0.0,1.0])
        ]
        rotation_list = [
            Rotation([0.0,0.0,1.0],2)
        ]


        improper_rotation_list = [
        ]

        label = 'C_{2h}'

        super().__init__(identity_list+inversion_list+reflection_list+rotation_list+improper_rotation_list, label)



class C2v(PointGroup):

    def __init__(self):
        identity_list = [
            Identity()
        ]

        inversion_list = [
        ]

        reflection_list = [
            Reflection([1.0,0.0,0.0]),
            Reflection([0.0,1.0,0.0])
        ]

        rotation_list = [
            Rotation([0.0,0.0,1.0],2)
        ]


        improper_rotation_list = [
        ]

        label = 'C_{2v}'

        super().__init__(identity_list+inversion_list+reflection_list+rotation_list+improper_rotation_list, label)


class D2(PointGroup):

    def __init__(self):
        identity_list = [
            Identity()
        ]

        inversion_list = [
        ]

        reflection_list = [
        ]

        rotation_list = [
            Rotation([1.0,0.0,0.0],2),
            Rotation([0.0,1.0,0.0],2),
            Rotation([0.0,0.0,1.0],2)
        ]


        improper_rotation_list = [
        ]

        label = 'D_{2}'

        super().__init__(identity_list+inversion_list+reflection_list+rotation_list+improper_rotation_list, label)

class C2(PointGroup):

    def __init__(self):
        identity_list = [
            Identity()
        ]

        inversion_list = [
        ]

        reflection_list = [
        ]

        rotation_list = [
            Rotation([0.0,0.0,1.0],2)
        ]


        improper_rotation_list = [
        ]

        label = 'C_{2}'

        super().__init__(identity_list+inversion_list+reflection_list+rotation_list+improper_rotation_list, label)


class Cs(PointGroup):

    def __init__(self):
        identity_list = [
            Identity()
        ]

        inversion_list = [
        ]

        reflection_list = [
            Reflection([0.0,0.0,1.0])
        ]

        rotation_list = [
        ]


        improper_rotation_list = [
        ]

        label = 'C_{s}'

        super().__init__(identity_list+inversion_list+reflection_list+rotation_list+improper_rotation_list, label)

class Ci(PointGroup):

    def __init__(self):
        identity_list = [
            Identity()
        ]

        inversion_list = [
            Inversion()
        ]

        reflection_list = [
        ]

        rotation_list = [
        ]


        improper_rotation_list = [
        ]

        label = 'C_{i}'

        super().__init__(identity_list+inversion_list+reflection_list+rotation_list+improper_rotation_list, label)

class C1(PointGroup):

    def __init__(self):
        identity_list = [
            Identity()
        ]

        inversion_list = [
        ]

        reflection_list = [
        ]

        rotation_list = [
        ]


        improper_rotation_list = [
        ]

        label = 'C_{1}'

        super().__init__(identity_list+inversion_list+reflection_list+rotation_list+improper_rotation_list, label)

