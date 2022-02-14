import numpy as np


class SZ:
    """non-spin-adapted spin label"""

    def __init__(self, n=0, twos=0, pg=0):
        self.n = n
        self.twos = twos
        self.pg = pg

    @property
    def is_fermion(self):
        return self.n % 2 == 1

    def __add__(self, other):
        return self.__class__(self.n + other.n, self.twos + other.twos, self.pg ^ other.pg)

    def __sub__(self, other):
        return self.__class__(self.n - other.n, self.twos - other.twos, self.pg ^ other.pg)

    def __neg__(self):
        return self.__class__(-self.n, -self.twos, self.pg)

    def __eq__(self, other):
        return self.n == other.n and self.twos == other.twos and self.pg == other.pg

    def __lt__(self, other):
        return (self.n, self.twos, self.pg) < (other.n, other.twos, other.pg)

    def __hash__(self):
        return hash((self.n, self.twos, self.pg))

    def __repr__(self):
        if self.twos % 2 == 1:
            return "< N=%d SZ=%d/2 PG=%d >" % (self.n, self.twos, self.pg)
        else:
            return "< N=%d SZ=%d PG=%d >" % (self.n, self.twos // 2, self.pg)

