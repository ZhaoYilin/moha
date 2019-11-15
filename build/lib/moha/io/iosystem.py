class IOSystem(object):
    def __init__(self):
        pass
    @classmethod
    def from_file(cls,geofile,basisfile):
        if geofile.endswith('.xyz'):
            from moha.io.iogeometry import load_xyz
            molecule = load_xyz(geofile)
        else:
            raise ValueError('Unknown file format')

            
        symbols = []
        coordinates = []
        for atom in molecule.atoms:
            symbols.append(atom.symbol)
            coordinates.append(atom.coordinate)

        if basisfile.endswith('.nwchem'):
            from moha.io.iobasis import load_nwchem
            basis_set = load_nwchem(symbols,coordinates,basisfile)
        else:
            raise ValueError('Unknown file format')

        return molecule,basis_set

