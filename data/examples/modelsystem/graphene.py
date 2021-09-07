from moha import *
from math import sqrt

a = 0.24595     # unit cell length
a_cc = 0.142    # carbon carbon distance

def graphene(a,a_cc):


    #define a 2D graphene lattice cell with vectors a1 a2 and a3
    cell = Cell(2,[[a,0.,0.],[a/2,a/2*sqrt(3),0.],[0.,0.,0.]])


    #add a site labeled 'A' at positon [0., -a_cc/2., 0.]
    cell.add_site(LatticeSite([0., -a_cc/2., 0.],'A'))
    #add a site labeled 'B' at positon [0., a_cc/2., 0.]
    cell.add_site(LatticeSite([0., a_cc/2., 0.],'B'))

    #add a bond from first site in cell [0,0,0] to second site in cell [0,0,0]
    cell.add_bond(LatticeBond([0,0,0],0,[0,0,0],1))
    #add a bond from second site in cell [0,0,0] to first site in cell [0,1,0]
    cell.add_bond(LatticeBond([0,0,0],1,[0,1,0],0))
    #add a bond from first site in cell [0,0,0] to second site in cell [1,-1,0]
    cell.add_bond(LatticeBond([0,0,0],0,[1,-1,0],1))

    #buld a lattice of shape [4,4,1] with cell we defined
    lattice = Lattice(cell,[4,4,1])
    return lattice

lattice = graphene(a,a_cc)
#plot the lattice we constructed
lattice.plot()
