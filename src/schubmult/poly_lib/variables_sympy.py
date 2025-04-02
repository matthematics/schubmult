from functools import cache

from sympy import Expr, Integer
from sympy.core.symbol import Str

import schubmult.poly_lib.variables as vsyme


class ISymbol(Expr):
    def __new__(cls, base, index):
        return ISymbol.__xnew_cached__(cls, base, index)

    @staticmethod
    @cache
    def __xnew_cached__(_class, base, index):
        return ISymbol.__xnew__(_class, base, index)

    @staticmethod
    def __xnew__(_class, base, index):
        return Expr.__new__(_class, base, Integer(index))

    def _sympystr(self, printer):
        return printer._print_Symbol(self)

    @property
    def base(self):
        return self.args[0]

    @property
    def name(self):
        return f"{self.args[0].label}_{self.args[1]}"

    @property
    def index(self):
        return self.args[1]

    @property
    def indices(self):
        return (self.args[1],)
    
    def __hash__(self):
        return hash(tuple(self.args))

    def __eq__(self, other):
        if not vsyme.is_indexed(other):
            return False
        return other.base == self.base and other.index == self.index
 
    def _symengine_(self):
        return vsyme.ISymbol(self.name, *self.args)
    