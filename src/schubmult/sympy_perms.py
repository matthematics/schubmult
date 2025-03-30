# from typing import NamedTuple

import sympy.combinatorics.permutations as spp

import schubmult.perm_lib as pl


class Permutation:
    def __init__(self, perm, sperm=None):
        p = tuple(pl.permtrim([*perm]))
        self._perm = p
        if sperm:
            self._sperm = sperm
        else:
            self._sperm = spp.Permutation._af_new([i-1 for i in p])

    def __getitem__(self, i):
        # print("yay")
        return self._perm[i]

    def __setitem__(self, i, v):
        raise NotImplementedError

    def __hash__(self):
        return hash(self._sperm)

    def __mul__(self, other):
        # print("yay")
        new_sperm =  other._sperm * self._sperm
        new_perm = pl.permtrim([new_sperm(i) + 1 for i in new_sperm.list()])
        if len(new_perm) != new_sperm.size:
            new_sperm = spp.Permutation._af_new([i-1 for i in new_perm])
        return Permutation(new_perm,new_sperm)

    def __iter__(self):
        # print("yay")
        return self._perm.__iter__()

    def __getslice__(self, i, j):
        # print("yay")
        return self._perm[i:j]

    def __str__(self):
        # print("yay")
        return str(self._perm)

    def __add__(self, other):
        # print("yay")
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*self._perm] + other
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    def __radd__(self, other):
        # print("yay")
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = other + [*self._perm]
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    def __eq__(self, other):
        # print("yay")
        if isinstance(other, Permutation):
            return other._perm == self._perm
        if isinstance(other, list):
            return [*self._perm] == other
        if isinstance(other, tuple):
            return self._perm == other
        return False

    def __len__(self):
        # print("yay")
        return len(self._perm)

    def __invert__(self):
        new_sperm = ~(self._sperm)
        new_perm = [new_sperm(i) + 1 for i in range(new_sperm.size)]
        return Permutation(new_perm, new_sperm)
    
    def __repr__(self):
        return self.__str__()

def inv(perm):
    return perm.inversions()

#def permtrim


