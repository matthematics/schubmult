# from typing import NamedTuple

from functools import cache, cached_property

import sympy.combinatorics.permutations as spp
from sympy import Basic

import schubmult.perm_lib as pl


class Permutation(Basic):

    def __new__(cls, perm, sperm=None):
        return Permutation.__xnew_cached__(cls, tuple(perm), sperm)

    @staticmethod
    @cache
    def __xnew_cached__(_class, perm, sperm):
        return Permutation.__xnew__(_class, perm, sperm)

    @staticmethod
    def __xnew__(_class, perm, sperm):
        obj = Basic.__new__(_class)
        if isinstance(perm, Permutation):
            # print("this is happening")
            obj._perm = perm._perm
            obj._sperm = perm._sperm
        else:
            p = tuple(pl.permtrim_list([*perm]))
            obj._perm = p
            if len(obj._perm)<2:
                obj._perm = (1,2)
            if sperm:
                obj._sperm = sperm
            else:
                obj._sperm = spp.Permutation._af_new([i - 1 for i in p])
        return obj


    @property
    def code(self):
        return self._sperm.inversion_vector()


    def inv(self):
        return self._sperm.inversions()

    def swap(self, i, j):
        import sys
        new_perm = [*self._perm]
        #print(f"{new_perm=}",file=sys.stderr)
        if i>j:
            i, j = j, i
        #print(f"OLD {new_perm=} {new_perm[i]=} {new_perm[j]=} {i=} {j=} FWIPO",file=sys.stderr)
        
        if j>=len(new_perm):
            new_perm += list(range(len(new_perm)+1,j+2))
            #print(f"bugs {new_perm=}", file=sys.stderr)
        new_perm[i], new_perm[j] = new_perm[j], new_perm[i]
        #print(f"NEW {new_perm=} {new_perm[i]=} {new_perm[j]=} FWoPO",file=sys.stderr)
        return Permutation(new_perm)


    def __getitem__(self, i):
        if isinstance(i, slice):
            return [self[ii] for ii in range(*i.indices(len(self._perm)))]
        if i >= len(self._perm):
            return i + 1
        return self._perm[i]

    def __setitem__(self, i, v):
        raise NotImplementedError

    def __hash__(self):
        return hash(self._perm)

    def __mul__(self, other):
        # print("yay")
        new_sperm = other._sperm * self._sperm
        new_perm = pl.permtrim_list([new_sperm(i) + 1 for i in range(new_sperm.size)])
        # print(f"{new_perm=}")
        if len(new_perm) != new_sperm.size:
            new_sperm = spp.Permutation._af_new([i - 1 for i in new_perm])
        return Permutation(new_perm, new_sperm)

    def __iter__(self):
        yield from self._perm.__iter__()

    def __getslice__(self, i, j):
        return self._perm[i:j]

    def __str__(self):
        # print("yay")
        return str(self._perm)

    def __add__(self, other):
        # print("yay")
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*self._perm, *other]
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    def __radd__(self, other):
        # print("yay")
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*other, *self._perm]
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


# def permtrim
