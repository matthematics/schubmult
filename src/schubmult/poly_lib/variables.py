# class generators with base
# symbols cls argument!

from functools import cache

import sympy
from symengine import Mul, symbols, sympify
from symengine.lib.symengine_wrapper import Symbol
from sympy.core.symbol import Str

import schubmult.poly_lib.variables_sympy as vsym


class ISymbol(Symbol):
    def __new__(cls, dummy_name, base, index):
        #print(f"frofulating bagel {name=} {base=} {index=}")
        return ISymbol.__xnew_cached__(cls, base, index)
    
    @staticmethod
    @cache
    def __xnew_cached__(_class, base, index):
        return ISymbol.__xnew__(_class, base, index)

    @staticmethod
    def __xnew__(_class, base, index):
        obj = Symbol.__new__(_class, f"{base.label}_{index}")
        obj._base = base
        obj._index = index
        return obj

    @property
    def base(self):
        return self._base

    @property
    def name(self):
        return f"{self._base.label}_{self._index}"

    @property
    def index(self):
        return self._index

    def __hash__(self):
        return hash((self._base,self._index))

    def __eq__(self, other):
        try:
            return self.name == sympify(other).name
        except Exception:
            return NotImplemented

    def _sympy_(self):
        return vsym.ISymbol(self._base, self._index)

    @property
    def args(self):
        return (self.name, self._base, self._index)

    @property
    def func(self):
        return self.__class__

class GeneratingSet(Str):
    def __new__(cls, name):
        return GeneratingSet.__xnew_cached__(cls, name)
    
    is_Atom = True

    @staticmethod
    @cache
    def __xnew_cached__(_class, name):
        return GeneratingSet.__xnew__(_class, name)

    @staticmethod
    def __xnew__(_class, name):
        obj = Str.__new__(_class, name)
        obj._symbols_arr = tuple([symbols(f"{name}_{i}",cls=ISymbol,base=obj,index=i) for i in range(100)]) 
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


class SMul(Mul,ISymbol):
    
    def __new__(cls, *args):
        obj = Mul.__new__(cls,*args)
        return obj
    
    def _sympy_(self):
        return vsym.SMul(*[sympy.sympify(arg) for arg in self.args])


def is_indexed(x):
    return hasattr(x, "index") and hasattr(x, "base")
