from moha import *

d = 1.0     # unit cell length

def cube(d):
    #define a 3D cubic lattice cell with vectors a1 a2 and a3
    cell = Cell(3,[d, 0., 0.],[0., d, 0.],[0., 0., d])

    #add a site labeled 'A' at positon [0.,0.,0.]
    cell.add_site(LatticeSite([0.,0.,0.],'A'))

    #add a bond from first site in cell [0,0,0] to second site in cell [0,0,0]
    cell.add_bond(LatticeBond([0,0,0],0,[1,0,0],0))
    #add a bond from second site in cell [0,0,0] to first site in cell [0,1,0]
    cell.add_bond(LatticeBond([0,0,0],0,[0,1,0],0))
    #add a bond from second site in cell [0,0,0] to first site in cell [0,1,0]
    cell.add_bond(LatticeBond([0,0,0],0,[0,0,1],0))

    #buld a lattice of shape [4,4,4] with cell we defined
    lattice = Lattice(cell,[4,4,4])
    return lattice

lattice = cube(d)
#plot the lattice we constructed
lattice.plot()
