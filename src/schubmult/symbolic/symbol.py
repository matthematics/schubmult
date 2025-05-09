from functools import cache

from sympy import Symbol

from .expr import SymengineExpr, SympyExpr


class SymengineSymbol(SymengineExpr):

    def __new__(cls, name):
        return SymengineSymbol.__xnew_cached__(cls, name)
    
    @staticmethod
    @cache
    def __xnew_cached__(_class, name):
        return SymengineSymbol.__xnew__(_class, name)
    
    @staticmethod
    def __xnew__(_class, name):
        return SymengineExpr.__new__(_class, name)

    def _sympystr(self, printer):
        return printer._print_Symbol(self)

    @cache
    def _sympy_(self):
        return SympySymbol(self)

    @property
    def name(self):
        return self.args[0]
    
    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if hasattr(other, "name"):
            return self.name == other.name
        return False

class SympySymbol(SympyExpr, Symbol):
    pass