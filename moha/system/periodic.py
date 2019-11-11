class Element(object):
    '''Represents an element from the periodic table.

       The following attributes are supported for all elements:

       number
            The atomic number.

       symbol
            A string with the symbol of the element.

       name
            The full element name.

       group
            The group of the element (not for actinides and lanthanides).

       period
            The row of the periodic system.
    '''

    def __init__(self, number=None, symbol=None, **kwargs):
        self.number = number
        self.symbol = symbol
        for name, value in kwargs.items():
            setattr(self, name, value)



class Periodic(object):
    '''A periodic table data structure.'''
    def __init__(self, elements):
        '''**Arguments:**

           elements
                A list of :class:`Element` instances.
        '''
        self.elements = elements
        self._lookup = {}
        for element in elements:
            self._lookup[element.number] = element
            self._lookup[element.symbol.lower()] = element

    def __getitem__(self, index):
        '''Get an element from the table based on a flexible index.

           **Argument:**

           index
                This can be either an integer atomic number, a string with the
                elemental symbol (any case), or a string with the atomic number.

           **Returns:** the corresponding :class:`Element` instance
        '''
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
    elements = []
    H = Element(1,'H',name ='hydrogen')
    He = Element(2,'He',name ='helium')
    Li = Element(3,'Li',name ='lithium')
    Be = Element(4,'Be',name ='Beryllium')
    B = Element(5,'B',name ='boron')
    C = Element(6,'C',name ='carbon')
    N = Element(7,'N',name ='nitrogen')
    O = Element(8,'O',name ='oxygen')
    F = Element(9,'F',name ='fluorine')
    Ne = Element(10,'Ne',name ='neon')
    Na = Element(11,'Na',name ='sodium')
    Mg = Element(12,'Mg',name ='magnesium')
    Al = Element(13,'Al',name ='aluminium')
    Si = Element(14,'Si',name ='silicon')
    P = Element(15,'P',name ='phosphorus')
    S = Element(16,'S',name ='sulfur')
    Cl = Element(17,'Cl',name ='chlorine')
    Ar = Element(18,'Ar',name ='argon')
    elements.append(H)
    elements.append(He)
    elements.append(Li)
    elements.append(Be)
    elements.append(B)
    elements.append(C)
    elements.append(N)
    elements.append(O)
    elements.append(F)
    elements.append(Ne)
    elements.append(Na)
    elements.append(Mg)
    elements.append(Al)
    elements.append(Si)
    elements.append(P)
    elements.append(S)
    elements.append(Cl)
    elements.append(Ar)
    return Periodic(elements)

