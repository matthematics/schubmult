# class generators with base
# symbols cls argument!

from functools import cache

import sympy
from symengine import Symbol, symbols, sympify
from sympy.core.symbol import Str


# variable registry
class GeneratingSet(Str):
    def __new__(cls, name):
        return GeneratingSet.__xnew_cached__(cls, name)

    _registry = {}

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
    if (v in GeneratingSet._registry) or (isinstance(v, sympy.Symbol) and sympify(v) in GeneratingSet._registry):
        return GeneratingSet._registry[v]
    if isinstance(sympify(v), Symbol):
        return sympify(v).name.split("_")[0], int(v.name.split("_")[1])
    raise ValueError(f"Unknown type: {type(v)}")
