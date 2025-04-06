# class generators with base
# symbols cls argument!

import re
from bisect import bisect_left
from functools import cache
from typing import ClassVar

from symengine import Symbol, symbols, sympify
from sympy import Basic
from sympy.core.symbol import Str


# variable registry
# TODO: ensure sympifies
# TODO: masked generating set
class GeneratingSet(Str):
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
        return GeneratingSet.__xnew__(_class, name)

    @staticmethod
    def __xnew__(_class, name):
        obj = Str.__new__(_class, name)
        obj._symbols_arr = tuple([symbols(f"{name}_{i}") for i in range(100)])
        obj._index_lookup = {obj._symbols_arr[i]: i for i in range(len(obj._symbols_arr))}
        return obj

    @property
    def label(self):
        return self.name

    # index of v in the genset
    def index(self, v):
        return self._index_lookup[v]

    def __repr__(self):
        return f"GeneratingSet('{self.label}')"

    def __str__(self):
        return self.name

    def _latex(self, printer):
        return printer._print_Str(self)

    def _sympystr(self, printer):
        return printer.doprint(self.name)

    def __getitem__(self, i):
        return self._symbols_arr[i]

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return isinstance(other, GeneratingSet) and self.label == other.label

class MaskedGeneratingSet(Basic):

    def __new__(cls, gset, index_mask):
        return MaskedGeneratingSet.__xnew_cached__(cls, gset, index_mask)

    @staticmethod
    @cache
    def __xnew_cached__(_class, gset, index_mask):
        return MaskedGeneratingSet.__xnew__(_class, gset, index_mask)

    @staticmethod
    def __xnew__(_class, gset, index_mask):
        obj = Basic.__new__(_class, gset, index_mask)
        # obj._symbols_arr = tuple([symbols(f"{name}_{i}") for i in range(100)])
        # obj._index_lookup = {obj._symbols_arr[i]: i for i in range(len(obj._symbols_arr))}
        mask_dict = {}
        cur_index = 0
        for i in range(len(gset._symbols_arr)):
            if bisect_left(index_mask, i) != i:
                mask_dict[cur_index] = i
                cur_index += 1
        obj._mask = mask_dict
        obj._index_lookup = {gset[index_mask[i]]: i for i in range(len(gset._symbols_arr) - len(index_mask))}
        return obj

    @property
    def base_genset(self):
        return self.args[0]

    def __getitem__(self, index):
        return self.base_genset[self._mask[index]]

    def index(self, v):
        return self._index_lookup(v)


def base_index(v):
    if isinstance(v, (list, tuple)):
        return base_index(v[0])[0], None
    if isinstance(v, GeneratingSet):
        return v.label, None
    try:
        v = sympify(v)
    except Exception:
        try:
            import sage  # noqa: F401
        except ImportError:
            return None, None
        from sage.rings.polynomial.multivariate_polynomial import MPolynomial  # type: ignore

        if isinstance(v, MPolynomial):
            m = GeneratingSet._sage_index_pattern.match(str(v))
            if m:
                return m.group(1), int(m.group(2))
    if v in GeneratingSet._registry:
        return GeneratingSet._registry[v]
    if isinstance(v, Symbol):
        try:
            m = GeneratingSet._index_pattern.match(v.name)
            if m:
                return m.group(1), int(m.group(2))
        except Exception:
            pass
        m = GeneratingSet._sage_index_pattern.match(str(v))
        if m:
            return m.group(1), int(m.group(2))
    return None, None
