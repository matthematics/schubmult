import math
from functools import cache, cached_property

import sympy.combinatorics.permutations as spp
from symengine import sympify
from sympy import Basic, Tuple

import schubmult.utils.logging as lg
from schubmult.utils.perm_utils import cyclic_sort, permtrim_list, trimcode

# schubmult.poly_lib.variables import GeneratingSet

logger = lg.get_logger(__name__)

zero = sympify(0)
n = 100

# TODO: permutations act


class Permutation(Basic):
    def __new__(cls, perm):
        return Permutation.__xnew_cached__(cls, tuple(perm))

    print_as_code = False

    @staticmethod
    @cache
    def __xnew_cached__(_class, perm):
        return Permutation.__xnew__(_class, perm)

    @staticmethod
    def __xnew__(_class, perm):
        p = tuple(permtrim_list([*perm]))
        s_perm = spp.Permutation._af_new([i - 1 for i in p])
        obj = Basic.__new__(_class, Tuple(*perm))
        obj._s_perm = s_perm
        obj._perm = p
        obj._hash_code = hash(p)
        cd = s_perm.inversion_vector()
        obj._unique_key = (len(p), sum([cd[i] * math.factorial(len(p) - 1 - i) for i in range(len(cd))]))
        return obj

    @classmethod
    def sorting_perm(cls, itera):
        L = [i + 1 for i in range(len(itera))]
        L.sort(key=lambda i: itera[i - 1])
        return Permutation(L)

    def _sympystr(self, printer):
        if Permutation.print_as_code:
            return printer.doprint(trimcode(self))
        return printer.doprint(self._perm)

    def __call__(self, i):
        """1-indexed"""
        return self[i - 1]

    def descents(self, zero_indexed=True):
        if zero_indexed:
            return self._s_perm.descents()
        return {i + 1 for i in self._s_perm.descents()}

    def get_cycles(self):
        return self.get_cycles_cached()

    @cache
    def get_cycles_cached(self):
        return [tuple(cyclic_sort([i + 1 for i in c])) for c in self._s_perm.cyclic_form]

    @property
    def code(self):
        return list(self.cached_code())

    @cache
    def cached_code(self):
        return self._s_perm.inversion_vector()

    @cached_property
    def inv(self):
        return self._s_perm.inversions()

    def swap(self, i, j):
        new_perm = [*self._perm]
        # print(f"SWAP {new_perm=}")
        if i > j:
            i, j = j, i
        if j >= len(new_perm):
            # print(f"SWAP {j}>={new_perm=}")
            new_perm += list(range(len(new_perm) + 1, j + 2))
            # print(f"SWAP extended {new_perm=}")
        new_perm[i], new_perm[j] = new_perm[j], new_perm[i]
        # print(f"SWAP iddle {new_perm=}")
        return Permutation(new_perm)

    def __getitem__(self, i):
        if isinstance(i, slice):
            return [self[ii] for ii in range(i.start if i.start is not None else 0, i.stop if i.stop is not None else len(self))]
        if i >= len(self._perm):
            return i + 1
        return self._perm[i]

    def __setitem__(self, i, v):
        raise NotImplementedError

    def __hash__(self):
        return self._hash_code

    def __mul__(self, other):
        new_sperm = other._s_perm * self._s_perm
        new_perm = permtrim_list([new_sperm.array_form[i] + 1 for i in range(new_sperm.size)])
        return Permutation(new_perm)

    def __iter__(self):
        yield from self._perm.__iter__()

    def __getslice__(self, i, j):
        return self._perm[i:j]

    def __str__(self):
        return str(self._perm)

    def __add__(self, other):
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*self._perm, *other]
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    def __radd__(self, other):
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*other, *self._perm]
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    def __eq__(self, other):
        if isinstance(other, Permutation):
            # print(f"{other._perm= } {self._perm=} {type(self._perm)=}")
            # return other._perm == self._perm
            return other._unique_key == self._unique_key
        if isinstance(other, list):
            # print(f"{[*self._perm]= } {other=}")
            return [*self._perm] == other
        if isinstance(other, tuple):
            # print(f"{self._perm=} {other=}")
            return self._perm == other
        return False

    def __len__(self):
        # print("REMOVE THIS")
        return max(len(self._perm), 2)

    def __invert__(self):
        new_sperm = ~(self._s_perm)
        new_perm = [new_sperm.array_form[i] + 1 for i in range(new_sperm.size)]
        return Permutation(new_perm)

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return tuple(self) < tuple(other)
