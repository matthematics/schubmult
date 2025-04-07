# class generators with base
# symbols cls argument!

import re
from bisect import bisect_left
from functools import cache
from typing import ClassVar

from symengine import Symbol, symbols, sympify
from sympy import Basic, Tuple
from sympy.core.symbol import Str


class GeneratingSet_base(Basic):
    def __new__(cls, *args):
        return Basic.__new__(cls, *args)

    def __getitem__(self, i):
        return NotImplemented

    def __len__(self):
        return NotImplemented

# variable registry
# TODO: ensure sympifies
# TODO: masked generating set
class GeneratingSet(GeneratingSet_base):
    def __new__(cls, name):
        return GeneratingSet.__xnew_cached__(cls, name)

    _registry: ClassVar = {}

    _index_pattern = re.compile("^([^_]+)_([0-9]+)$")
    _sage_index_pattern = re.compile("^([^0-9]+)([0-9]+)$")
    # is_Atom = True
    # TODO: masked generating set
    @staticmethod
    @cache
    def __xnew_cached__(_class, name):
        return GeneratingSet.__xnew__(_class, Str(str(name)))

    @staticmethod
    def __xnew__(_class, name):
        obj = GeneratingSet_base.__new__(_class, name)
        obj._symbols_arr = tuple([symbols(f"{name}_{i}") for i in range(100)])
        obj._index_lookup = {obj._symbols_arr[i]: i for i in range(len(obj._symbols_arr))}
        return obj

    #@property

    @property
    def label(self):
        return str(self.args[0])

    # index of v in the genset
    def index(self, v):
        return self._index_lookup[v]

    def __repr__(self):
        return f"GeneratingSet('{self.label}')"

    def __str__(self):
        return self.name

    def _latex(self, printer):
        return printer.doprint(self.label)

    def _sympystr(self, printer):
        return printer.doprint(self.label)

    def __getitem__(self, i):
        return self._symbols_arr[i]

    def __len__(self):
        return len(self._symbols_arr)

    def __hash__(self):
        return hash(self.label)

    def __iter__(self):
        yield from [self[i] for i in range(len(self))]

    def __eq__(self, other):
        return isinstance(other, GeneratingSet) and self.label == other.label


class MaskedGeneratingSet(GeneratingSet_base):
    def __new__(cls, gset, index_mask):
        return MaskedGeneratingSet.__xnew_cached__(cls, gset, tuple(sorted(index_mask)))

    @staticmethod
    @cache
    def __xnew_cached__(_class, gset, index_mask):
        return MaskedGeneratingSet.__xnew__(_class, gset, index_mask)

    @staticmethod
    def __xnew__(_class, gset, index_mask):
        obj = GeneratingSet_base.__new__(_class, gset, Tuple(*index_mask))
        # obj._symbols_arr = tuple([symbols(f"{name}_{i}") for i in range(100)])
        # obj._index_lookup = {obj._symbols_arr[i]: i for i in range(len(obj._symbols_arr))}
        mask_dict = {}
        mask_dict[0] = 0
        cur_index = 1
        for i in range(1,len(gset._symbols_arr)):
            index = bisect_left(index_mask, i)
            if index>=len(index_mask) or index_mask[index] != i:
                mask_dict[cur_index] = i
                cur_index += 1
        # print(f"{index_mask=} {mask_dict=}")
        obj._mask = mask_dict
        obj._index_lookup = {gset[i]: mask_dict[i] for i in range(len(gset) - len(index_mask))}
        obj._label = gset.label
        return obj

    @property
    def label(self):
        return str(self._label)

    def set_label(self, label):
        self._label = label

    @property
    def index_mask(self):
        return tuple(self.args[1])

    def complement(self):
        return MaskedGeneratingSet(self.base_genset, [i for i in range(1,len(self.base_genset)) if i not in set(self.index_mask)])

    @property
    def base_genset(self):
        return self.args[0]

    def __getitem__(self, index):
        return self.base_genset[self._mask[index]]

    def __iter__(self):
        yield from [self[i] for i in range(len(self))]

    def index(self, v):
        return self._index_lookup(v)

    def __len__(self):
        return len(self.base_genset) - len(self.index_mask)

# def base_index(v):
#     if isinstance(v, (list, tuple)):
#         return base_index(v[0])[0], None
#     if isinstance(v, GeneratingSet):
#         return v.label, None
#     try:
#         v = sympify(v)
#     except Exception:
#         try:
#             import sage  # noqa: F401
#         except ImportError:
#             return None, None
#         from sage.rings.polynomial.multivariate_polynomial import MPolynomial  # type: ignore

#         if isinstance(v, MPolynomial):
#             m = GeneratingSet._sage_index_pattern.match(str(v))
#             if m:
#                 return m.group(1), int(m.group(2))
#     if v in GeneratingSet._registry:
#         return GeneratingSet._registry[v]
#     if isinstance(v, Symbol):
#         try:
#             m = GeneratingSet._index_pattern.match(v.name)
#             if m:
#                 return m.group(1), int(m.group(2))
#         except Exception:
#             pass
#         m = GeneratingSet._sage_index_pattern.match(str(v))
#         if m:
#             return m.group(1), int(m.group(2))
#     return None, None
