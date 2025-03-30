# from typing import NamedTuple

from functools import cached_property

import sympy.combinatorics.permutations as spp

import schubmult.perm_lib as pl


class Permutation:
    def __init__(self, perm, sperm=None):
        if isinstance(perm, Permutation):
            # print("this is happening")
            self._perm = perm._perm
            self._sperm = perm._sperm
        else:
            p = tuple(pl.permtrim_list([*perm]))
            self._perm = p
            if sperm:
                self._sperm = sperm
            else:
                self._sperm = spp.Permutation._af_new([i - 1 for i in p])

    @property
    def code(self):
        return self._sperm.inversion_vector()

    def __getitem__(self, i):
        if isinstance(i, slice):
            return [self[ii] for ii in range(*i.indices(len(self._perm)))]
        if i >= len(self._perm):
            return i + 1
        return self._perm[i]

    def __setitem__(self, i, v):
        raise NotImplementedError

    def __hash__(self):
        return hash(self._sperm)

    def __mul__(self, other):
        # print("yay")
        new_sperm = other._sperm * self._sperm
        new_perm = pl.permtrim_list([new_sperm(i) + 1 for i in range(new_sperm.size)])
        # print(f"{new_perm=}")
        if len(new_perm) != new_sperm.size:
            new_sperm = spp.Permutation._af_new([i - 1 for i in new_perm])
        return Permutation(new_perm, new_sperm)

    def __iter__(self):
        return self._perm.__iter__()

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
