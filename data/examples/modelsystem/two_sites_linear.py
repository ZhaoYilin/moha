from moha import *

d = 1.0     # unit cell length
d_ab = 0.2  # distance between site A and B

def linear(d,d_ab):
    #define a 1D linear lattice cell with vectors a1 a2 and a3
    cell = Cell(1,[d, 0., 0.],[0., 0., 0.],[0., 0., 0.])
    #add a site labeled 'A' at positon [0., 0., 0.]
    cell.add_site(LatticeSite([0., 0., 0.], 'A'))
    #add a site labeled 'B' at positon [d_ab, 0., 0.]
    cell.add_site(LatticeSite([d_ab, 0., 0.],'B'))
    #add a bond from first site in cell [0,0,0] to second site in cell [0,0,0]
    cell.add_bond(LatticeBond([0,0,0],0,[0,0,0],1))
    #add a bond from second site in cell [0,0,0] to first site in cell [1,0,0]
    cell.add_bond(LatticeBond([0,0,0],1,[1,0,0],0))
    #buld a lattice of shape [4,1,1] with cell we defined
    lattice = Lattice(cell,[4,1,1])
    return lattice

lattice = linear(d,d_ab)
#plot the lattice we constructed
lattice.plot() 
