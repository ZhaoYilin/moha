from moha import *
import sys

class Molden(object):
    """
    """
    def __init__(self,filename,mol,orbs,wf):
        """
        """
        self.filename = filename
        self.mol = mol
        self.orbs = orbs
        self.wf = wf

    def dump_molden(self):
        """write molden input file
        
        **Arguments:**

        filename
            File name of molden input file

        data
        """
        with open(self.filename+'.molden','w') as f:
            self.write_header(f)
            self.write_atoms(f)
            self.write_GTO(f)
            self.write_MO(f)

    def write_header(self,handler):
        handler.write("[Molden Format]\n")
        handler.write('\n')

    def write_atoms(self,handler):
        handler.write("[Atoms] AU\n")
        for i,atom in enumerate(self.mol.atoms):
            handler.write("%s \t %d \t  %d \t %f \t %f \t %f\n" 
            % (atom.symbol,i+1,atom.number,atom.coordinate[0],atom.coordinate[1],atom.coordinate[2]))

    def write_GTO(self,handler):
        handler.write("[GTO]\n")

        atom_set_index = self.atom_index_filter()
        orbital_set_index = self.orbital_index_filter()
        for i,index in enumerate(orbital_set_index):
            if index in atom_set_index[1:]:
                handler.write("\n")
            if index in atom_set_index:
                handler.write("%d \t 0\n"%(atom_set_index.index(index)+1))
            self.write_atom_orbital_header(handler,self.orbs.bases[index])
            for j in range(len(self.orbs.bases[index].exps)):
                handler.write("\t\t\t %f \t %f \n"%(self.orbs.bases[index].exps[j],self.orbs.bases[index].coefs[j]))

    def write_MO(self,handler):
        handler.write("[MO]\n")
        for i in range(self.wf.nspatial):
            self.write_MO_header(handler,i)
            self.write_MO_coefficient(handler,i)

    
    def atom_index_filter(self):
        """filte the orbital
        """
        list = []
        for i,basis in enumerate(self.orbs.bases):
            list.append(basis.atom_index)

        inds, unq = [],[]
        seen = set()
        for i,ele in enumerate(list):
            if ele not in seen:
                inds.append(i)
                unq.append(ele)
            seen.add(ele) 
        return inds
    
    def orbital_index_filter(self):
        """filte the orbital
        """
        list = []
        for i,basis in enumerate(self.orbs.bases):
            tmp = []
            tmp.append(basis.atom_index)
            tmp.append(basis.n_number)
            tmp.append(sum(basis.shell))
            list.append(tuple(tmp))

        inds, unq = [],[]
        seen = set()
        for i,ele in enumerate(list):
            if ele not in seen:
                inds.append(i)
                unq.append(ele)
            seen.add(ele) 
        return inds
    
    def write_atom_orbital_header(self,handler,basis):
        if sum(basis.shell) == 0:
            handler.write("s \t %d \t 1.00\n"%(len(basis.exps)))
        elif sum(basis.shell) == 1:
            handler.write("p \t %d \t 1.00\n"%(len(basis.exps)))
        elif sum(basis.shell) == 2:
            handler.write("d \t %d \t 1.00\n"%(len(basis.exps)))
        elif sum(basis.shell) == 3:
            handler.write("f \t %d \t 1.00\n"%(len(basis.exps)))
        elif sum(basis.shell) == 4:
            handler.write("g \t %d \t 1.00\n"%(len(basis.exps)))
        else:
            print('error')

    def write_MO_header(self,handler,i):
        occ_alpha = [1]*self.wf.occ['alpha']+[0]*(self.wf.nspatial-self.wf.occ['alpha'])
        occ_beta = [1]*self.wf.occ['beta']+[0]*(self.wf.nspatial-self.wf.occ['beta'])
        handler.write("Sym=%s\n"%('A'))
        handler.write("Ene=%f\n"%(self.wf.orbital_energies[i]))
        handler.write("Spin=%s\n"%('Alpha'))
        handler.write("Occup=%f\n"%(occ_alpha[i]+occ_beta[i]))

    def write_MO_coefficient(self,handler,i):
        for j in range(self.wf.nspatial):
            handler.write("{} \t {} \n".format(j+1,self.wf.coefficients[j,i]))
