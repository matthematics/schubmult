# class generators with base
# symbols cls argument!

import re
from functools import cache

from symengine import Symbol, symbols, sympify
from sympy.core.symbol import Str


# variable registry
class GeneratingSet(Str):
    def __new__(cls, name):
        return GeneratingSet.__xnew_cached__(cls, name)

    _registry = {}

    _index_pattern = re.compile("^([^_]+)_([0-9]+)$")
    is_Atom = True

    @staticmethod
    @cache
    def __xnew_cached__(_class, name):
        return GeneratingSet.__xnew__(_class, name)

    @staticmethod
    def __xnew__(_class, name):
        obj = Str.__new__(_class, name)
        obj._symbols_arr = tuple([symbols(f"{name}_{i}") for i in range(100)])
        for i in range(100):
            GeneratingSet._registry[obj._symbols_arr[i]] = (name, i)
        return obj

    @property
    def label(self):
        return self.name

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



def base_index(v):
    v = sympify(v)
    if v in GeneratingSet._registry:
        return GeneratingSet._registry[v]
    if isinstance(v, Symbol):
        m = GeneratingSet._index_pattern.match(v.name)
        if m:
            return m.group(1), int(m.group(2))
    return NotImplemented, NotImplemented
