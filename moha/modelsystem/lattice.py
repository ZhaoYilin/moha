'''Lattice for Model Hamiltonian'''
import copy
import numpy as np
import matplotlib.pyplot as plt

site_label_dir = {'0':'A','1':'B','2':'C','3':'D','4':'E'}
site_color_dir = {'0':'b','1':'r','2':'g'}
bond_label_dir = {'0':'a','1':'b','2':'c','3':'d','4':'e'}
bond_color_dir = {'0':'b','1':'r','2':'g'}

class LatticeSite(object):
    """Child class of LatticeSite"""

    def __init__(self,coordinate,label=None,color=None,markersize=80,marker='o'):
        """ Initialize of lattice  site instance

        Parameters
        ----------
        coordinate : list
            coordinate of the site

        label : str
            label of the site

        color : str
            color the the site in matplotlib

        marker : str
            type of marker in matplotlib

        markersize : int
            size of the marker in matplotlib

        """
        self.coordinate = np.array(coordinate)
        self.label = label
        self.color = color
        self.marker = marker
        self.markersize = markersize


class LatticeBond(object):
    """Bond of Lattice"""

    def __init__(self,cell1,site1,cell2,site2,color=None,label=None,linewidth=None,linetype=None):
        """ Initialize of lattice bond instance

        Parameters
        ----------
        cell1 : list
            integer coordinate of first cell

        site1 : int
            index of first site

        cell2 : list
            integer coordinate of second cell

        site2 : int
            index of second site

        label : str
            label of the bond

        color : str
            color of the line in matplotlib

        linewidth : int
            integer to represent the width of line in matplotlib

        linetype : str
            str to represent the type of line in matplotlib

        """
        self.cell1 = np.array(cell1)
        self.site1 = site1
        self.cell2 = np.array(cell2)
        self.site2 = site2
        coordinate1 = None
        coordinate2 = None
        self.color = color
        self.label = label
        self.linewidth = linewidth
        self.linetype = linetype


class Cell(object):
    """Primitive cell of the lattice"""

    def __init__(self,dimension,a1,a2,a3):
        """ Initialize of cell instance

        Parameters
        ----------
        dimension : int
            dimension of the cell, which can be 1, 2 or 3

        a1 : list
            first component of primitive vector 

        a2 : list
            second component of primitive vector 

        a3 : list
            third component of primitive vector 

        """
        self.dimension = dimension
        self.a1 = np.array(a1)
        self.a2 = np.array(a2)
        self.a3 = np.array(a3)
        self.sites = []
        self.bonds = []
        self.Nsites = 0
        self.Nbonds = 0 

    def add_site(self,site):
        """Add a lattice site instance to the cell

        Parameters
        ----------

        site : LatticeXSite 
            lattice site class

        """
        if site.label==None:
            site.label = site_label_dir[str(self.Nsites)]
        if site.color==None:
            site.color = site_color_dir[str(self.Nsites)]
        self.sites.append(site)
        self.Nsites +=1

    def add_bond(self,bond):
        """Add a lattice bond instance to the cell

        Parameters
        ----------

        bond : LatticeBond 
            lattice bond class

        """
        if bond.label==None:
            bond.label = bond_label_dir[str(self.Nbonds)]
        if bond.color==None:
            bond.color = 'black'
        self.bonds.append(bond)
        self.Nbonds +=1


class Lattice(object):
    """Class for lattice"""

    def __init__(self, cell, shape):
        """Initialize of lattice instance

        Parameters
        ----------

        cell : Cell 
            Cell class

        shape : list
            a list of three integer

        dimension : int
            dimension of the lattice, which can be one two or three

        Nsites : int
            number of sites in the lattice

        sites : numpy array
            numpy array for element of site object

        bonds : list
            list of bond object

        """
        self.cell = cell
        self.shape = shape
        self.dimension = cell.dimension
        self.Nsites = np.prod(shape) * self.cell.Nsites
        self.sites = np.zeros(self.shape+[self.cell.Nsites],dtype='object')
        self.bonds = []
        self.build_sites()
        self.build_bonds()

    def build_sites(self):
        """Add all sites to the lattice."""
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                for k in range(self.shape[2]):
                    for s,site in enumerate(self.cell.sites):
                        newsite = copy.deepcopy(site)
                        coordinate = self.cell.a1*i+\
                        self.cell.a2*j+\
                        self.cell.a3*k
                        newsite.coordinate += coordinate
                        self.sites[i,j,k,s] = newsite

    def build_bonds(self):
        """Add all bonds to the lattice."""
        shape_prime = np.array([self.shape[0]-1,self.shape[1]-1,self.shape[2]-1])
        zeros = np.array([0,0,0])
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                for k in range(self.shape[2]):
                    for b,bond in enumerate(self.cell.bonds):
                        newbond = copy.deepcopy(bond)
                        newbond.cell1 += [i,j,k]
                        newbond.cell2 += [i,j,k]
                        #ToDo make a function to shorten those lines
                        if np.prod(newbond.cell1 <= shape_prime) and np.prod(newbond.cell2<=shape_prime) and np.prod(zeros <=newbond.cell1) and np.prod(zeros <= newbond.cell2):
                            newbond.coordinate1 = self.sites[newbond.cell1[0],newbond.cell1[1],newbond.cell1[2],newbond.site1].coordinate
                            newbond.coordinate2 = self.sites[newbond.cell2[0],newbond.cell2[1],newbond.cell2[2],newbond.site2].coordinate
                            self.bonds.append(newbond)

    def plotsite(self):
        """plot sites with matplotlib."""
        if self.dimension==1 or self.dimension==2:
            for site in self.sites.flatten():
                plt.scatter(site.coordinate[0],site.coordinate[1],marker=site.marker,s=site.markersize,c=site.color)
        elif self.dimension==3:
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure()
            ax = Axes3D(fig)
            for site in self.sites.flatten():
                ax.scatter(site.coordinate[0],site.coordinate[1],site.coordinate[2],marker=site.marker,s=site.markersize,c=site.color)

    def plotbond(self):
        """plot bonds with matplotlib."""
        if self.dimension==1 or self.dimension==2:
            for bond in self.bonds:
                x = [bond.coordinate1[0],bond.coordinate2[0]]
                y = [bond.coordinate1[1],bond.coordinate2[1]]
                plt.plot(x,y,c=bond.color)
        elif self.dimension==3:
            for bond in self.bonds:
                x = [bond.coordinate1[0],bond.coordinate2[0]]
                y = [bond.coordinate1[1],bond.coordinate2[1]]
                z = [bond.coordinate1[2],bond.coordinate2[2]]

                plt.plot(x,y,z,c=bond.color)

    def plot(self):
        """plot lattice and show it with matplotlib."""
        self.plotsite()
        self.plotbond()
        plt.show()
