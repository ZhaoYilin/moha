# Periodic table used in MoHa
# Default unit using in the peridoic table
# mass: amu
# radius: pm

class Element(object):
    """Represents an element from the periodic table.
               
    Attributes
    ----------
    number : int
        The atomic number.

    symbol : str
        A string with the symbol of the element.

    name : str
        The full element name.

    group : int
        The group of the element (not for actinides and lanthanides).
    """    

    def __init__(self, number=None, symbol=None, **kwargs):
        self.number = number
        self.symbol = symbol
        for name, value in kwargs.items():
            setattr(self, name, value)



class Periodic(object):
    """A periodic table data structure.
    """
    def __init__(self, elements):
        """Initialize the instance

        Parameters
        ----------
        elements
            A list of :class:`Element` instances.
        """
        self.elements = elements
        self._lookup = {}
        for element in elements:
            self._lookup[element.number] = element
            self._lookup[element.symbol.lower()] = element

    def __getitem__(self, index):
        """Get an element from the table based on a flexible index.

        Parameters
        ----------
        index
            This can be either an integer atomic number, a string with the
            elemental symbol (any case), or a string with the atomic number.

        Returns
        -------
        result : 
            The corresponding :class:`Element` instance.
        """
        result = self._lookup.get(index)
        if result is None and isinstance(index, str):
            index = index.strip()
            result = self._lookup.get(index.lower())
            if result is None and index.isdigit():
                result = self._lookup.get(int(index))
                if result is None:
                    raise KeyError('Could not find element %s.' % index)
        return result

def load_periodic():
    """Load the periodic table data.
       
    Unit set
    --------
    radius_calculated: pm
        

    Returns
    -------
    periodic_table:
        An instance of :class:`Periodic`.
        
    """
    elements = []
    #First row
    H = Element(1,'H',name ='Hydrogen',period=1,group=1,mass=1.007975,radius_calculated=53,radius_covalent=37,radius_van_der_waals=120)
    He = Element(2,'He',name ='Helium',period=1,group=18,mass=4.002602,radius_calculated=31,radius_covalent=32,radius_van_der_waals=140)
    #Second row
    Li = Element(3,'Li',name ='Lithium',period=2,group=1,mass=6.9675,radius_calculated=167,radius_covalent=134,radius_van_der_waals=182)
    Be = Element(4,'Be',name ='Beryllium',period=2,group=2,mass=9.0121831,radius_calculated=112,radius_covalent=90,radius_van_der_waals=None)
    B = Element(5,'B',name ='Boron',period=2,group=13,mass=10.8135,radius_calculated=87,radius_covalent=82,radius_van_der_waals=None)
    C = Element(6,'C',name ='Carbon',period=2,group=14,mass=12.0106,radius_calculated=67,radius_covalent=77,radius_van_der_waals=170)
    N = Element(7,'N',name ='Nitrogen',period=2,group=15,mass=14.006855,radius_calculated=56,radius_covalent=75,radius_van_der_waals=155)
    O = Element(8,'O',name ='Oxygen',period=2,group=16,mass=15.9994,radius_calculated=48,radius_covalent=73,radius_van_der_waals=152)
    F = Element(9,'F',name ='Fluorine',period=2,group=17,mass=18.998403163,radius_calculated=42,radius_covalent=71,radius_van_der_waals=147)
    Ne = Element(10,'Ne',name ='Neon',period=2,group=18,mass=20.1797,radius_calculated=38,radius_covalent=69,radius_van_der_waals=154)
    #Third row
    Na = Element(11,'Na',name ='Sodium',period=3,group=1,mass=22.98976928,radius_calculated=190,radius_covalent=154,radius_van_der_waals=227)
    Mg = Element(12,'Mg',name ='Magnesium',period=3,group=2,mass=24.3055,radius_calculated=145,radius_covalent=130,radius_van_der_waals=173)
    Al = Element(13,'Al',name ='Aluminium',period=3,group=13,mass=26.9815385,radius_calculated=118,radius_covalent=118,radius_van_der_waals=None)
    Si = Element(14,'Si',name ='Silicon',period=3,group=14,mass=28.085,radius_calculated=111,radius_covalent=111,radius_van_der_waals=210)
    P = Element(15,'P',name ='Phosphorus',period=3,group=15,mass=30.973761998,radius_calculated=98,radius_covalent=106,radius_van_der_waals=180)
    S = Element(16,'S',name ='Sulfur',period=3,group=16,mass=32.0675,radius_calculated=88,radius_covalent=102,radius_van_der_waals=180)
    Cl = Element(17,'Cl',name ='Chlorine',period=3,group=17,mass=35.4515,radius_calculated=79,radius_covalent=99,radius_van_der_waals=175)
    Ar = Element(18,'Ar',name ='Argon',period=3,group=18,mass=39.948,radius_calculated=71,radius_covalent=97,radius_van_der_waals=188)
    #Fourth row
    K = Element(19,'K',name ='Potassium',period=4,group=1,mass=39.0983,radius_calculated=243,radius_covalent=196,radius_van_der_waals=275)
    Ca = Element(20,'Ca',name ='Calcium',period=4,group=2,mass=40.078,radius_calculated=194,radius_covalent=174,radius_van_der_waals=None)
    Sc = Element(21,'Sc',name ='Scandium',period=4,group=3,mass=44.955908,radius_calculated=184,radius_covalent=144,radius_van_der_waals=None)
    Ti = Element(22,'Ti',name ='Titanium',period=4,group=4,mass=47.867,radius_calculated=176,radius_covalent=136,radius_van_der_waals=None)
    V = Element(23,'V',name ='Vanadium',period=4,group=5,mass=50.9415,radius_calculated=171,radius_covalent=125,radius_van_der_waals=None)
    Cr = Element(24,'Cr',name ='Chromium',period=4,group=6,mass=51.9961,radius_calculated=166,radius_covalent=127,radius_van_der_waals=None)
    Mn = Element(25,'Mn',name ='Manganese',period=4,group=7,mass=54.938044,radius_calculated=161,radius_covalent=139,radius_van_der_waals=None)
    Fe = Element(26,'Fe',name ='Iron',period=4,group=8,mass=55.845,radius_calculated=156,radius_covalent=125,radius_van_der_waals=None)
    Co = Element(27,'Co',name ='Cobalt',period=4,group=9,mass=58.933194,radius_calculated=152,radius_covalent=126,radius_van_der_waals=None)
    Ni = Element(28,'Ni',name ='Nickel',period=4,group=10,mass=58.6934,radius_calculated=149,radius_covalent=121,radius_van_der_waals=163)
    Cu = Element(29,'Cu',name ='Copper',period=4,group=11,mass=63.546,radius_calculated=145,radius_covalent=138,radius_van_der_waals=140)
    Zn = Element(30,'Zn',name ='Zinc',period=4,group=12,mass=65.38,radius_calculated=142,radius_covalent=131,radius_van_der_waals=139)
    Ga = Element(31,'Ga',name ='Gallium',period=4,group=13,mass=69.723,radius_calculated=136,radius_covalent=126,radius_van_der_waals=187)
    Ge = Element(32,'Ge',name ='Germanium',period=4,group=14,mass=72.63,radius_calculated=125,radius_covalent=122,radius_van_der_waals=None)
    As = Element(33,'As',name ='Arsenic',period=4,group=15,mass=74.921595,radius_calculated=114,radius_covalent=119,radius_van_der_waals=185)
    Se = Element(34,'Se',name ='Selenium',period=4,group=16,mass=78.971,radius_calculated=103,radius_covalent=116,radius_van_der_waals=190)
    Br = Element(35,'Br',name ='Bromine',period=4,group=17,mass=79.904,radius_calculated=94,radius_covalent=114,radius_van_der_waals=185)
    Kr = Element(36,'Kr',name ='Krypton',period=4,group=18,mass=83.798,radius_calculated=88,radius_covalent=110,radius_van_der_waals=202)
    #Fifth row
    Rb = Element(37,'Rb',name ='Rubidium',period=5,group=1,mass=85.4678,radius_calculated=265,radius_covalent=211,radius_van_der_waals=None)
    Sr = Element(38,'Sr',name ='Strontium',period=5,group=2,mass=87.62,radius_calculated=219,radius_covalent=192,radius_van_der_waals=None)
    Y = Element(39,'Y',name ='Yttrium',period=5,group=3,mass=88.90584,radius_calculated=212,radius_covalent=162,radius_van_der_waals=None)
    Zr = Element(40,'Zr',name ='Zirconium',period=5,group=4,mass=91.224,radius_calculated=206,radius_covalent=148,radius_van_der_waals=None)
    Nb = Element(41,'Nb',name ='Niobium',period=5,group=5,mass=92.90637,radius_calculated=198,radius_covalent=137,radius_van_der_waals=None)
    Mo = Element(42,'Mo',name ='Molybdenum',period=5,group=6,mass=95.95,radius_calculated=190,radius_covalent=145,radius_van_der_waals=None)
    Tc = Element(43,'Tc',name ='Technetium',period=5,group=7,mass=96.90636,radius_calculated=183,radius_covalent=156,radius_van_der_waals=None)
    Ru = Element(44,'Ru',name ='Ruthenium',period=5,group=8,mass=101.07,radius_calculated=178,radius_covalent=126,radius_van_der_waals=None)
    Rh = Element(45,'Rh',name ='Rhodium',period=5,group=9,mass=102.9055,radius_calculated=173,radius_covalent=135,radius_van_der_waals=None)
    Pd = Element(46,'Pd',name ='Palladium',period=5,group=10,mass=106.42,radius_calculated=169,radius_covalent=131,radius_van_der_waals=163)
    Ag = Element(47,'Ag',name ='Silver',period=5,group=11,mass=107.8682,radius_calculated=165,radius_covalent=153,radius_van_der_waals=172)
    Cd = Element(48,'Cd',name ='Cadmium',period=5,group=12,mass=112.414,radius_calculated=161,radius_covalent=148,radius_van_der_waals=158)
    In = Element(49,'In',name ='Indium',period=5,group=13,mass=114.818,radius_calculated=156,radius_covalent=144,radius_van_der_waals=193)
    Sn = Element(50,'Sn',name ='Tin',period=5,group=14,mass=118.71,radius_calculated=145,radius_covalent=141,radius_van_der_waals=217)
    Sb = Element(51,'Sb',name ='Antimony',period=5,group=15,mass=121.76,radius_calculated=133,radius_covalent=138,radius_van_der_waals=None)
    Te = Element(52,'Te',name ='Tellurium',period=5,group=16,mass=127.6,radius_calculated=123,radius_covalent=135,radius_van_der_waals=206)
    I = Element(53,'I',name ='Iodine',period=5,group=17,mass=126.90447,radius_calculated=115,radius_covalent=133,radius_van_der_waals=198)
    Xe = Element(54,'Xe',name ='Xenon',period=5,group=18,mass=131.293,radius_calculated=108,radius_covalent=130,radius_van_der_waals=216)
    #Sixth row
    Cs = Element(55,'Cs',name ='Cesium',period=6,group=1,mass=132.90545196,radius_calculated=298,radius_covalent=225,radius_van_der_waals=None)
    Ba = Element(56,'Ba',name ='Barium',period=6,group=2,mass=137.327,radius_calculated=253,radius_covalent=198,radius_van_der_waals=None)
    La = Element(57,'La',name ='Lanthanum',period=6,group=3,mass=138.90547,radius_calculated=None,radius_covalent=169,radius_van_der_waals=None)
    Ce = Element(58,'Ce',name ='Cerium',period=6,group=3,mass=140.116,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Pr = Element(59,'Pr',name ='Praseodymium',period=6,group=3,mass=140.90766,radius_calculated=247,radius_covalent=None,radius_van_der_waals=None)
    Nd = Element(60,'Nd',name ='Neodymium',period=6,group=3,mass=144.242,radius_calculated=206,radius_covalent=None,radius_van_der_waals=None)
    Pm = Element(61,'Pm',name ='Promethium',period=6,group=3,mass=144.91276,radius_calculated=205,radius_covalent=None,radius_van_der_waals=None)
    Sm = Element(62,'Sm',name ='Samarium',period=6,group=3,mass=150.36,radius_calculated=238,radius_covalent=None,radius_van_der_waals=None)
    Eu = Element(63,'Eu',name ='Europium',period=6,group=3,mass=151.964,radius_calculated=231,radius_covalent=None,radius_van_der_waals=None)
    Gd = Element(64,'Gd',name ='Gadolinium',period=6,group=3,mass=157.25,radius_calculated=233,radius_covalent=None,radius_van_der_waals=None)
    Tb = Element(65,'Tb',name ='Terbium',period=6,group=3,mass=158.92535,radius_calculated=225,radius_covalent=None,radius_van_der_waals=None)
    Dy = Element(66,'Dy',name ='Dysprosium',period=6,group=3,mass=162.5,radius_calculated=228,radius_covalent=None,radius_van_der_waals=None)
    Ho = Element(67,'Ho',name ='Holmium',period=6,group=3,mass=164.93033,radius_calculated=226,radius_covalent=None,radius_van_der_waals=None)
    Er = Element(68,'Er',name ='Erbium',period=6,group=3,mass=167.259,radius_calculated=226,radius_covalent=None,radius_van_der_waals=None)
    Tm = Element(69,'Tm',name ='Thulium',period=6,group=3,mass=168.93422,radius_calculated=222,radius_covalent=None,radius_van_der_waals=None)
    Yb = Element(70,'Yb',name ='Ytterbium',period=6,group=3,mass=173.054,radius_calculated=222,radius_covalent=None,radius_van_der_waals=None)
    Lu = Element(71,'Lu',name ='Lutetium',period=6,group=3,mass=174.9668,radius_calculated=217,radius_covalent=160,radius_van_der_waals=None)
    Hf = Element(72,'Hf',name ='Hafnium',period=6,group=4,mass=178.49,radius_calculated=208,radius_covalent=150,radius_van_der_waals=None)
    Ta = Element(73,'Ta',name ='Tantalum',period=6,group=5,mass=180.94788,radius_calculated=200,radius_covalent=138,radius_van_der_waals=None)
    W = Element(74,'W',name ='Tungsten',period=6,group=6,mass=183.84,radius_calculated=193,radius_covalent=146,radius_van_der_waals=None)
    Re = Element(75,'Re',name ='Rhenium',period=6,group=7,mass=186.207,radius_calculated=188,radius_covalent=159,radius_van_der_waals=None)
    Os = Element(76,'Os',name ='Osmium',period=6,group=8,mass=190.23,radius_calculated=185,radius_covalent=128,radius_van_der_waals=None)
    Ir = Element(77,'Ir',name ='Iridium',period=6,group=9,mass=192.217,radius_calculated=180,radius_covalent=137,radius_van_der_waals=None)
    Pt = Element(78,'Pt',name ='Platinum',period=6,group=10,mass=195.084,radius_calculated=177,radius_covalent=128,radius_van_der_waals=175)
    Au = Element(79,'Au',name ='Gold',period=6,group=11,mass=196.966569,radius_calculated=174,radius_covalent=144,radius_van_der_waals=166)
    Hg = Element(80,'Hg',name ='Mercury',period=6,group=12,mass=200.592,radius_calculated=171,radius_covalent=149,radius_van_der_waals=155)
    Tl = Element(81,'Tl',name ='Thallium',period=6,group=13,mass=204.3835,radius_calculated=156,radius_covalent=148,radius_van_der_waals=196)
    Pb = Element(82,'Pb',name ='Lead',period=6,group=14,mass=207.2,radius_calculated=154,radius_covalent=147,radius_van_der_waals=202)
    Bi = Element(83,'Bi',name ='Bismuth',period=6,group=15,mass=208.9804,radius_calculated=143,radius_covalent=146,radius_van_der_waals=None)
    Po = Element(84,'Po',name ='Polonium',period=6,group=16,mass=208.98243,radius_calculated=135,radius_covalent=None,radius_van_der_waals=None)
    At = Element(85,'At',name ='Astatine',period=6,group=17,mass=209.98715,radius_calculated=127,radius_covalent=None,radius_van_der_waals=None)
    Rn = Element(86,'Rn',name ='Radon',period=6,group=18,mass=222.01758,radius_calculated=120,radius_covalent=145,radius_van_der_waals=None)
    #Seventh row
    Fr = Element(87,'Fr',name ='Francium',period=7,group=1,mass=223.01973,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Ra = Element(88,'Ra',name ='Radium',period=7,group=2,mass=226.02541,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Ac = Element(89,'Ac',name ='Actinium',period=7,group=3,mass=227.02775,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Th = Element(90,'Th',name ='Thorium',period=7,group=3,mass=232.038,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Pa = Element(91,'Pa',name ='Protactinium',period=7,group=3,mass=231.03588,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    U = Element(92,'U',name ='Uranium',period=7,group=3,mass=238.0289,radius_calculated=None,radius_covalent=None,radius_van_der_waals=186)
    Np = Element(93,'Np',name ='Neptunium',period=7,group=3,mass=237.048172,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Pu = Element(94,'Pu',name ='Plutonium',period=7,group=3,mass=244.06420,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Am = Element(95,'Am',name ='Americium',period=7,group=3,mass=243.061380,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Cm = Element(96,'Gm',name ='Curium',period=7,group=3,mass=247.07035,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Bk = Element(97,'Bk',name ='Berkelium',period=7,group=3,mass=247.07031,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Cf = Element(98,'Cf',name ='Californium',period=7,group=3,mass=251.07959,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Es = Element(99,'Es',name ='Einsteinium',period=7,group=3,mass=252.0830,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Fm = Element(100,'Fm',name ='Fermium',period=7,group=3,mass=257.09511,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Md = Element(101,'Md',name ='Mendelevium',period=7,group=3,mass=258.09843,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    No = Element(102,'No',name ='Nobelium',period=7,group=3,mass=259.10100,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Lr = Element(103,'Lr',name ='Lawrencium',period=7,group=3,mass=266.120,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Rf = Element(104,'Rf',name ='Rutherfordium',period=7,group=4,mass=267.122,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Db = Element(105,'Db',name ='Dubnium',period=7,group=5,mass=268.126,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Sg = Element(106,'Sg',name ='Seaborgium',period=7,group=6,mass=269.128,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Bh = Element(107,'Bh',name ='Bohrium',period=7,group=7,mass=270.133,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Hs = Element(108,'Hs',name ='Hassium',period=7,group=8,mass=269.1336,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Mt = Element(109,'Mt',name ='Meitnerium',period=7,group=9,mass=277.154,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Ds = Element(110,'Ds',name ='Darmstadtium',period=7,group=10,mass=282.166,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Rg = Element(111,'Rg',name ='Roentgenium',period=7,group=11,mass=282.169,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Cn = Element(112,'Cn',name ='Copernicium',period=7,group=12,mass=286.179,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Nh = Element(113,'Nh',name ='Nihonium',period=7,group=13,mass=286.182,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Fl = Element(114,'Fl',name ='Flerovium',period=7,group=14,mass=290.192,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Mc = Element(115,'Mc',name ='Moscovium',period=7,group=15,mass=290.196,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Lv = Element(116,'Lv',name ='Livermorium',period=7,group=16,mass=293.205,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Ts = Element(117,'Ts',name ='Tennessine',period=7,group=17,mass=294.211,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    Og = Element(118,'Og',name ='Oganesson',period=7,group=18,mass=295.216,radius_calculated=None,radius_covalent=None,radius_van_der_waals=None)
    #First row
    elements.append(H)
    elements.append(He)
    #Second row
    elements.append(Li)
    elements.append(Be)
    elements.append(B)
    elements.append(C)
    elements.append(N)
    elements.append(O)
    elements.append(F)
    elements.append(Ne)
    #Third row
    elements.append(Na)
    elements.append(Mg)
    elements.append(Al)
    elements.append(Si)
    elements.append(P)
    elements.append(S)
    elements.append(Cl)
    elements.append(Ar)
    #Fourth row
    elements.append(K)
    elements.append(Ca)
    elements.append(Sc)
    elements.append(Ti)
    elements.append(V)
    elements.append(Cr)
    elements.append(Mn)
    elements.append(Fe)
    elements.append(Co)
    elements.append(Ni)
    elements.append(Cu)
    elements.append(Zn)
    elements.append(Ga)
    elements.append(Ge)
    elements.append(As)
    elements.append(Se)
    elements.append(Br)
    elements.append(Kr)
    #Fifth row
    elements.append(Rb)
    elements.append(Sr)
    elements.append(Y)
    elements.append(Zr)
    elements.append(Nb)
    elements.append(Mo)
    elements.append(Tc)
    elements.append(Ru)
    elements.append(Rh)
    elements.append(Pd)
    elements.append(Ag)
    elements.append(Cd)
    elements.append(In)
    elements.append(Sn)
    elements.append(Sb)
    elements.append(Te)
    elements.append(I)
    elements.append(Xe)
    #Sixth row
    elements.append(Cs)
    elements.append(Ba)
    elements.append(La)
    elements.append(Ce)
    elements.append(Pr)
    elements.append(Nd)
    elements.append(Pm)
    elements.append(Sm)
    elements.append(Eu)
    elements.append(Gd)
    elements.append(Tb)
    elements.append(Dy)
    elements.append(Ho)
    elements.append(Er)
    elements.append(Tm)
    elements.append(Yb)
    elements.append(Lu)
    elements.append(Hf)
    elements.append(Ta)
    elements.append(W)
    elements.append(Re)
    elements.append(Os)
    elements.append(Ir)
    elements.append(Pt)
    elements.append(Au)
    elements.append(Hg)
    elements.append(Tl)
    elements.append(Pb)
    elements.append(Bi)
    elements.append(Po)
    elements.append(At)
    elements.append(Rn)
    #Seventh row
    elements.append(Fr)
    elements.append(Ra)
    elements.append(Ac)
    elements.append(Th)
    elements.append(Pa)
    elements.append(U)
    elements.append(Np)
    elements.append(Pu)
    elements.append(Am)
    elements.append(Cm)
    elements.append(Bk)
    elements.append(Cf)
    elements.append(Es)
    elements.append(Fm)
    elements.append(Md)
    elements.append(No)
    elements.append(Lr)
    elements.append(Rf)
    elements.append(Db)
    elements.append(Sg)
    elements.append(Bh)
    elements.append(Hs)
    elements.append(Mt)
    elements.append(Ds)
    elements.append(Rg)
    elements.append(Cn)
    elements.append(Nh)
    elements.append(Fl)
    elements.append(Mc)
    elements.append(Lv)
    elements.append(Ts)
    elements.append(Og)
    
    return Periodic(elements)

