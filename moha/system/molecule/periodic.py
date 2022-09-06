# Periodic table used in MoHa
# Default unit using in the peridoic table
# mass: amu
# radius: pm

__all__ = ['table']

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

    def __init__(self, **kwargs):
        for name, value in kwargs.items():
            setattr(self, name, value)


#Zeroth row
X = Element(
number=119,
symbol='X',
name ='Dummy')

Bq = Element(
number=120,
symbol='Bq',
name ='Ghost')

#First row
H = Element(
number=1,
symbol='H',
name ='Hydrogen',
period=1,
group=1,
mass=1.007975,
radius_calculated=53,
radius_covalent=37,
radius_van_der_waals=120)

He = Element(
number=2,
symbol='He',
name ='Helium',
period=1,
group=18,
mass=4.002602,
radius_calculated=31,
radius_covalent=32,
radius_van_der_waals=140)

#Second row
Li = Element(
number=3,
symbol='Li',
name ='Lithium',
period=2,
group=1,
mass=6.9675,
radius_calculated=167,
radius_covalent=134,
radius_van_der_waals=182)

Be = Element(
number=4,
symbol='Be',
name ='Beryllium',
period=2,
group=2,
mass=9.0121831,
radius_calculated=112,
radius_covalent=90,
radius_van_der_waals=None)

B = Element(
number=5,
symbol='B',
name ='Boron',
period=2,
group=13,
mass=10.8135,
radius_calculated=87,
radius_covalent=82,
radius_van_der_waals=None)

C = Element(
number=6,
symbol='C',
name ='Carbon',
period=2,
group=14,
mass=12.0106,
radius_calculated=67,
radius_covalent=77,
radius_van_der_waals=170)

N = Element(
number=7,
symbol='N',
name ='Nitrogen',
period=2,
group=15,
mass=14.006855,
radius_calculated=56,
radius_covalent=75,
radius_van_der_waals=155)

O = Element(
number=8,
symbol='O',
name ='Oxygen',
period=2,
group=16,
mass=15.9994,
radius_calculated=48,
radius_covalent=73,
radius_van_der_waals=152)

F = Element(
number=9,
symbol='F',
name ='Fluorine',
period=2,
group=17,
mass=18.998403163,
radius_calculated=42,
radius_covalent=71,
radius_van_der_waals=147)

Ne = Element(
number=10,
symbol='Ne',
name ='Neon',
period=2,
group=18,
mass=20.1797,
radius_calculated=38,
radius_covalent=69,
radius_van_der_waals=154)

#Third row
Na = Element(
number=11,
symbol='Na',
name ='Sodium',
period=3,
group=1,
mass=22.98976928,
radius_calculated=190,
radius_covalent=154,
radius_van_der_waals=227)

Mg = Element(
number=12,
symbol='Mg',
name ='Magnesium',
period=3,
group=2,
mass=24.3055,
radius_calculated=145,
radius_covalent=130,
radius_van_der_waals=173)

Al = Element(
number=13,
symbol='Al',
name ='Aluminium',
period=3,
group=13,
mass=26.9815385,
radius_calculated=118,
radius_covalent=118,
radius_van_der_waals=None)

Si = Element(
number=14,
symbol='Si',
name ='Silicon',
period=3,
group=14,
mass=28.085,
radius_calculated=111,
radius_covalent=111,
radius_van_der_waals=210)

P = Element(
number=15,
symbol='P',
name ='Phosphorus',
period=3,
group=15,
mass=30.973761998,
radius_calculated=98,
radius_covalent=106,
radius_van_der_waals=180)

S = Element(
number=16,
symbol='S',
name ='Sulfur',
period=3,
group=16,
mass=32.0675,
radius_calculated=88,
radius_covalent=102,
radius_van_der_waals=180)

Cl = Element(
number=17,
symbol='Cl',
name ='Chlorine',
period=3,
group=17,
mass=35.4515,
radius_calculated=79,
radius_covalent=99,
radius_van_der_waals=175)

Ar = Element(
number=18,
symbol='Ar',
name ='Argon',
period=3,
group=18,
mass=39.948,
radius_calculated=71,
radius_covalent=97,
radius_van_der_waals=188)

#Fourth row
K = Element(
number=19,
symbol='K',
name ='Potassium',
period=4,
group=1,
mass=39.0983,
radius_calculated=243,
radius_covalent=196,
radius_van_der_waals=275)

Ca = Element(
number=20,
symbol='Ca',
name ='Calcium',
period=4,
group=2,
mass=40.078,
radius_calculated=194,
radius_covalent=174,
radius_van_der_waals=None)

Sc = Element(
number=21,
symbol='Sc',
name ='Scandium',
period=4,
group=3,
mass=44.955908,
radius_calculated=184,
radius_covalent=144,
radius_van_der_waals=None)

Ti = Element(
number=22,
symbol='Ti',
name ='Titanium',
period=4,
group=4,
mass=47.867,
radius_calculated=176,
radius_covalent=136,
radius_van_der_waals=None)

V = Element(
number=23,
symbol='V',
name ='Vanadium',
period=4,
group=5,
mass=50.9415,
radius_calculated=171,
radius_covalent=125,
radius_van_der_waals=None)

Cr = Element(
number=24,
symbol='Cr',
name ='Chromium',
period=4,
group=6,
mass=51.9961,
radius_calculated=166,
radius_covalent=127,
radius_van_der_waals=None)

Mn = Element(
number=25,
symbol='Mn',
name ='Manganese',
period=4,
group=7,
mass=54.938044,
radius_calculated=161,
radius_covalent=139,
radius_van_der_waals=None)

Fe = Element(
number=26,
symbol='Fe',
name ='Iron',
period=4,
group=8,
mass=55.845,
radius_calculated=156,
radius_covalent=125,
radius_van_der_waals=None)

Co = Element(
number=27,
symbol='Co',
name ='Cobalt',
period=4,
group=9,
mass=58.933194,
radius_calculated=152,
radius_covalent=126,
radius_van_der_waals=None)

Ni = Element(
number=28,
symbol='Ni',
name ='Nickel',
period=4,
group=10,
mass=58.6934,
radius_calculated=149,
radius_covalent=121,
radius_van_der_waals=163)

Cu = Element(
number=29,
symbol='Cu',
name ='Copper',
period=4,
group=11,
mass=63.546,
radius_calculated=145,
radius_covalent=138,
radius_van_der_waals=140)

Zn = Element(
number=30,
symbol='Zn',
name ='Zinc',
period=4,
group=12,
mass=65.38,
radius_calculated=142,
radius_covalent=131,
radius_van_der_waals=139)

Ga = Element(
number=31,
symbol='Ga',
name ='Gallium',
period=4,
group=13,
mass=69.723,
radius_calculated=136,
radius_covalent=126,
radius_van_der_waals=187)

Ge = Element(
number=32,
symbol='Ge',
name ='Germanium',
period=4,
group=14,
mass=72.63,
radius_calculated=125,
radius_covalent=122,
radius_van_der_waals=None)

As = Element(
number=33,
symbol='As',
name ='Arsenic',
period=4,
group=15,
mass=74.921595,
radius_calculated=114,
radius_covalent=119,
radius_van_der_waals=185)

Se = Element(
number=34,
symbol='Se',
name ='Selenium',
period=4,
group=16,
mass=78.971,
radius_calculated=103,
radius_covalent=116,
radius_van_der_waals=190)

Br = Element(
number=35,
symbol='Br',
name ='Bromine',
period=4,
group=17,
mass=79.904,
radius_calculated=94,
radius_covalent=114,
radius_van_der_waals=185)

Kr = Element(
number=36,
symbol='Kr',
name ='Krypton',
period=4,
group=18,
mass=83.798,
radius_calculated=88,
radius_covalent=110,
radius_van_der_waals=202)

#Fifth row
Rb = Element(
number=37,
symbol='Rb',
name ='Rubidium',
period=5,
group=1,
mass=85.4678,
radius_calculated=265,
radius_covalent=211,
radius_van_der_waals=None)

Sr = Element(
number=38,
symbol='Sr',
name ='Strontium',
period=5,
group=2,
mass=87.62,
radius_calculated=219,
radius_covalent=192,
radius_van_der_waals=None)

Y = Element(
number=39,
symbol='Y',
name ='Yttrium',
period=5,
group=3,
mass=88.90584,
radius_calculated=212,
radius_covalent=162,
radius_van_der_waals=None)

Zr = Element(
number=40,
symbol='Zr',
name ='Zirconium',
period=5,
group=4,
mass=91.224,
radius_calculated=206,
radius_covalent=148,
radius_van_der_waals=None)

Nb = Element(
number=41,
symbol='Nb',
name ='Niobium',
period=5,
group=5,
mass=92.90637,
radius_calculated=198,
radius_covalent=137,
radius_van_der_waals=None)

Mo = Element(
number=42,
symbol='Mo',
name ='Molybdenum',
period=5,
group=6,
mass=95.95,
radius_calculated=190,
radius_covalent=145,
radius_van_der_waals=None)

Tc = Element(
number=43,
symbol='Tc',
name ='Technetium',
period=5,
group=7,
mass=96.90636,
radius_calculated=183,
radius_covalent=156,
radius_van_der_waals=None)

Ru = Element(
number=44,
symbol='Ru',
name ='Ruthenium',
period=5,
group=8,
mass=101.07,
radius_calculated=178,
radius_covalent=126,
radius_van_der_waals=None)

Rh = Element(
number=45,
symbol='Rh',
name ='Rhodium',
period=5,
group=9,
mass=102.9055,
radius_calculated=173,
radius_covalent=135,
radius_van_der_waals=None)

Pd = Element(
number=46,
symbol='Pd',
name ='Palladium',
period=5,
group=10,
mass=106.42,
radius_calculated=169,
radius_covalent=131,
radius_van_der_waals=163)

Ag = Element(
number=47,
symbol='Ag',
name ='Silver',
period=5,
group=11,
mass=107.8682,
radius_calculated=165,
radius_covalent=153,
radius_van_der_waals=172)

Cd = Element(
number=48,
symbol='Cd',
name ='Cadmium',
period=5,
group=12,
mass=112.414,
radius_calculated=161,
radius_covalent=148,
radius_van_der_waals=158)

In = Element(
number=49,
symbol='In',
name ='Indium',
period=5,
group=13,
mass=114.818,
radius_calculated=156,
radius_covalent=144,
radius_van_der_waals=193)

Sn = Element(
number=50,
symbol='Sn',
name ='Tin',
period=5,
group=14,
mass=118.71,
radius_calculated=145,
radius_covalent=141,
radius_van_der_waals=217)

Sb = Element(
number=51,
symbol='Sb',
name ='Antimony',
period=5,
group=15,
mass=121.76,
radius_calculated=133,
radius_covalent=138,
radius_van_der_waals=None)

Te = Element(
number=52,
symbol='Te',
name ='Tellurium',
period=5,
group=16,
mass=127.6,
radius_calculated=123,
radius_covalent=135,
radius_van_der_waals=206)

I = Element(
number=53,
symbol='I',
name ='Iodine',
period=5,
group=17,
mass=126.90447,
radius_calculated=115,
radius_covalent=133,
radius_van_der_waals=198)

Xe = Element(
number=54,
symbol='Xe',
name ='Xenon',
period=5,
group=18,
mass=131.293,
radius_calculated=108,
radius_covalent=130,
radius_van_der_waals=216)

#Sixth row
Cs = Element(
number=55,
symbol='Cs',
name ='Cesium',
period=6,
group=1,
mass=132.90545196,
radius_calculated=298,
radius_covalent=225,
radius_van_der_waals=None)

Ba = Element(
number=56,
symbol='Ba',
name ='Barium',
period=6,
group=2,
mass=137.327,
radius_calculated=253,
radius_covalent=198,
radius_van_der_waals=None)

La = Element(
number=57,
symbol='La',
name ='Lanthanum',
period=6,
group=3,
mass=138.90547,
radius_calculated=None,
radius_covalent=169,
radius_van_der_waals=None)

Ce = Element(
number=58,
symbol='Ce',
name ='Cerium',
period=6,
group=3,
mass=140.116,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Pr = Element(
number=59,
symbol='Pr',
name ='Praseodymium',
period=6,
group=3,
mass=140.90766,
radius_calculated=247,
radius_covalent=None,
radius_van_der_waals=None)

Nd = Element(
number=60,
symbol='Nd',
name ='Neodymium',
period=6,
group=3,
mass=144.242,
radius_calculated=206,
radius_covalent=None,
radius_van_der_waals=None)

Pm = Element(
number=61,
symbol='Pm',
name ='Promethium',
period=6,
group=3,
mass=144.91276,
radius_calculated=205,
radius_covalent=None,
radius_van_der_waals=None)

Sm = Element(
number=62,
symbol='Sm',
name ='Samarium',
period=6,
group=3,
mass=150.36,
radius_calculated=238,
radius_covalent=None,
radius_van_der_waals=None)

Eu = Element(
number=63,
symbol='Eu',
name ='Europium',
period=6,
group=3,
mass=151.964,
radius_calculated=231,
radius_covalent=None,
radius_van_der_waals=None)

Gd = Element(
number=64,
symbol='Gd',
name ='Gadolinium',
period=6,
group=3,
mass=157.25,
radius_calculated=233,
radius_covalent=None,
radius_van_der_waals=None)

Tb = Element(
number=65,
symbol='Tb',
name ='Terbium',
period=6,
group=3,
mass=158.92535,
radius_calculated=225,
radius_covalent=None,
radius_van_der_waals=None)

Dy = Element(
number=66,
symbol='Dy',
name ='Dysprosium',
period=6,
group=3,
mass=162.5,
radius_calculated=228,
radius_covalent=None,
radius_van_der_waals=None)

Ho = Element(
number=67,
symbol='Ho',
name ='Holmium',
period=6,
group=3,
mass=164.93033,
radius_calculated=226,
radius_covalent=None,
radius_van_der_waals=None)

Er = Element(
number=68,
symbol='Er',
name ='Erbium',
period=6,
group=3,
mass=167.259,
radius_calculated=226,
radius_covalent=None,
radius_van_der_waals=None)

Tm = Element(
number=69,
symbol='Tm',
name ='Thulium',
period=6,
group=3,
mass=168.93422,
radius_calculated=222,
radius_covalent=None,
radius_van_der_waals=None)

Yb = Element(
number=70,
symbol='Yb',
name ='Ytterbium',
period=6,
group=3,
mass=173.054,
radius_calculated=222,
radius_covalent=None,
radius_van_der_waals=None)

Lu = Element(
number=71,
symbol='Lu',
name ='Lutetium',
period=6,
group=3,
mass=174.9668,
radius_calculated=217,
radius_covalent=160,
radius_van_der_waals=None)

Hf = Element(
number=72,
symbol='Hf',
name ='Hafnium',
period=6,
group=4,
mass=178.49,
radius_calculated=208,
radius_covalent=150,
radius_van_der_waals=None)

Ta = Element(
number=73,
symbol='Ta',
name ='Tantalum',
period=6,
group=5,
mass=180.94788,
radius_calculated=200,
radius_covalent=138,
radius_van_der_waals=None)

W = Element(
number=74,
symbol='W',
name ='Tungsten',
period=6,
group=6,
mass=183.84,
radius_calculated=193,
radius_covalent=146,
radius_van_der_waals=None)

Re = Element(
number=75,
symbol='Re',
name ='Rhenium',
period=6,
group=7,
mass=186.207,
radius_calculated=188,
radius_covalent=159,
radius_van_der_waals=None)

Os = Element(
number=76,
symbol='Os',
name ='Osmium',
period=6,
group=8,
mass=190.23,
radius_calculated=185,
radius_covalent=128,
radius_van_der_waals=None)

Ir = Element(
number=77,
symbol='Ir',
name ='Iridium',
period=6,
group=9,
mass=192.217,
radius_calculated=180,
radius_covalent=137,
radius_van_der_waals=None)

Pt = Element(
number=78,
symbol='Pt',
name ='Platinum',
period=6,
group=10,
mass=195.084,
radius_calculated=177,
radius_covalent=128,
radius_van_der_waals=175)

Au = Element(
number=79,
symbol='Au',
name ='Gold',
period=6,
group=11,
mass=196.966569,
radius_calculated=174,
radius_covalent=144,
radius_van_der_waals=166)

Hg = Element(
number=80,
symbol='Hg',
name ='Mercury',
period=6,
group=12,
mass=200.592,
radius_calculated=171,
radius_covalent=149,
radius_van_der_waals=155)

Tl = Element(
number=81,
symbol='Tl',
name ='Thallium',
period=6,
group=13,
mass=204.3835,
radius_calculated=156,
radius_covalent=148,
radius_van_der_waals=196)

Pb = Element(
number=82,
symbol='Pb',
name ='Lead',
period=6,
group=14,
mass=207.2,
radius_calculated=154,
radius_covalent=147,
radius_van_der_waals=202)

Bi = Element(
number=83,
symbol='Bi',
name ='Bismuth',
period=6,
group=15,
mass=208.9804,
radius_calculated=143,
radius_covalent=146,
radius_van_der_waals=None)

Po = Element(
number=84,
symbol='Po',
name ='Polonium',
period=6,
group=16,
mass=208.98243,
radius_calculated=135,
radius_covalent=None,
radius_van_der_waals=None)

At = Element(
number=85,
symbol='At',
name ='Astatine',
period=6,
group=17,
mass=209.98715,
radius_calculated=127,
radius_covalent=None,
radius_van_der_waals=None)

Rn = Element(
number=86,
symbol='Rn',
name ='Radon',
period=6,
group=18,
mass=222.01758,
radius_calculated=120,
radius_covalent=145,
radius_van_der_waals=None)

#Seventh row
Fr = Element(
number=87,
symbol='Fr',
name ='Francium',
period=7,
group=1,
mass=223.01973,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Ra = Element(
number=88,
symbol='Ra',
name ='Radium',
period=7,
group=2,
mass=226.02541,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Ac = Element(
number=89,
symbol='Ac',
name ='Actinium',
period=7,
group=3,
mass=227.02775,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Th = Element(
number=90,
symbol='Th',
name ='Thorium',
period=7,
group=3,
mass=232.038,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Pa = Element(
number=91,
symbol='Pa',
name ='Protactinium',
period=7,
group=3,
mass=231.03588,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

U = Element(
number=92,
symbol='U',
name ='Uranium',
period=7,
group=3,
mass=238.0289,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=186)

Np = Element(
number=93,
symbol='Np',
name ='Neptunium',
period=7,
group=3,
mass=237.048172,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Pu = Element(
number=94,
symbol='Pu',
name ='Plutonium',
period=7,
group=3,
mass=244.06420,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Am = Element(
number=95,
symbol='Am',
name ='Americium',
period=7,
group=3,
mass=243.061380,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Cm = Element(
number=96,
symbol='Gm',
name ='Curium',
period=7,
group=3,
mass=247.07035,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Bk = Element(
number=97,
symbol='Bk',
name ='Berkelium',
period=7,
group=3,
mass=247.07031,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Cf = Element(
number=98,
symbol='Cf',
name ='Californium',
period=7,
group=3,
mass=251.07959,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Es = Element(
number=99,
symbol='Es',
name ='Einsteinium',
period=7,
group=3,
mass=252.0830,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Fm = Element(
number=100,
symbol='Fm',
name ='Fermium',
period=7,
group=3,
mass=257.09511,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Md = Element(
number=101,
symbol='Md',
name ='Mendelevium',
period=7,
group=3,
mass=258.09843,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

No = Element(
number=102,
symbol='No',
name ='Nobelium',
period=7,
group=3,
mass=259.10100,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Lr = Element(
number=103,
symbol='Lr',
name ='Lawrencium',
period=7,
group=3,
mass=266.120,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Rf = Element(
number=104,
symbol='Rf',
name ='Rutherfordium',
period=7,
group=4,
mass=267.122,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Db = Element(
number=105,
symbol='Db',
name ='Dubnium',
period=7,
group=5,
mass=268.126,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Sg = Element(
number=106,
symbol='Sg',
name ='Seaborgium',
period=7,
group=6,
mass=269.128,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Bh = Element(
number=107,
symbol='Bh',
name ='Bohrium',
period=7,
group=7,
mass=270.133,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Hs = Element(
number=108,
symbol='Hs',
name ='Hassium',
period=7,
group=8,
mass=269.1336,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Mt = Element(
number=109,
symbol='Mt',
name ='Meitnerium',
period=7,
group=9,
mass=277.154,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Ds = Element(
number=110,
symbol='Ds',
name ='Darmstadtium',
period=7,
group=10,
mass=282.166,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Rg = Element(
number=111,
symbol='Rg',
name ='Roentgenium',
period=7,
group=11,
mass=282.169,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Cn = Element(
number=112,
symbol='Cn',
name ='Copernicium',
period=7,
group=12,
mass=286.179,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Nh = Element(
number=113,
symbol='Nh',
name ='Nihonium',
period=7,
group=13,
mass=286.182,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Fl = Element(
number=114,
symbol='Fl',
name ='Flerovium',
period=7,
group=14,
mass=290.192,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Mc = Element(
number=115,
symbol='Mc',
name ='Moscovium',
period=7,
group=15,
mass=290.196,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Lv = Element(
number=116,
symbol='Lv',
name ='Livermorium',
period=7,
group=16,
mass=293.205,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Ts = Element(
number=117,
symbol='Ts',
name ='Tennessine',
period=7,
group=17,
mass=294.211,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

Og = Element(
number=118,
symbol='Og',
name ='Oganesson',
period=7,
group=18,
mass=295.216,
radius_calculated=None,
radius_covalent=None,
radius_van_der_waals=None)

table = {'':None,'h':H,'he':He,'li':Li,'be':Be,'b':B,'c':C,'n':N,'o':O,'f':F,'ne':Ne,'na':Na,'mg':Mg,'al':Al,'si':Si,'p':P,'s':S,'cl':Cl,'ar':Ar,'k':K,'ca':Ca,'sc':Sc,'ti':Ti,'v':V,'cr':Cr,'mc':Mn,'fe':Fe,'co':Co,'ni':Ni,'cu':Cu,'zn':Zn,'ga':Ga,'ge':Ge,'as':As,'se':Se,'br':Br,'kr':Kr,'rb':Rb,'sr':Sr,'y':Y,'zr':Zr,'nb':Nb,'mo':Mo,'tc':Tc,'ru':Ru,'rh':Rh,'pd':Pd,'ag':Ag,'cd':Cd,'in':In,'sn':Sn,'sb':Sb,'te':Te,'i':I,'xe':Xe,'cs':Cs,'ba':Ba,'la':La,'ce':Ce,'pr':Pr,'nd':Nd,'pm':Pm,'sm':Sm,'eu':Eu,'gd':Gd,'tb':Tb,'dy':Dy,'ho':Ho,'er':Er,'tm':Tm,'yb':Yb,'lu':Lu,'hf':Hf,'ta':Ta,'w':W,'re':Re,'os':Os,'ir':Ir,'pt':Pt,'au':Au,'hg':Hg,'tl':Tl,'pb':Pb,'bi':Bi,'po':Po,'at':At,'rn':Rn,'fr':Fr,'ra':Ra,'ac':Ac,'th':Th,'pa':Pa,'u':U,'np':Np,'pt':Pu,'am':Am,'cm':Cm,'bk':Bk,'cf':Cf,'es':Es,'fm':Fm,'md':Md,'no':No,'lr':Lr,'rf':Rf,'db':Db,'sg':Sg,'bh':Bh,'hs':Hs,'mt':Mt,'ds':Ds,'rg':Rg,'cn':Cn,'nh':Nh,'fl':Fl,'mc':Mc,'lv':Lv,'ts':Ts,'og':Og,'x':X,'bq':Bq}
