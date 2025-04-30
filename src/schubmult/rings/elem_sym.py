from functools import cache

from sympy import Integer, Tuple
from sympy.core.expr import Expr


class ElemSym(Expr):
    def __new__(cls, p, k, var1, var2):
        return ElemSym.__xnew_cached__(cls, p, k, tuple(var1), tuple(var2))

    @staticmethod
    @cache
    def __xnew_cached__(_class, p, k, var1, var2):
        return ElemSym.__xnew__(_class, p, k, var1, var2)

    @staticmethod
    def __xnew__(_class, p, k, var1, var2):
        obj = Expr.__new__(_class, Integer(p), Integer(k), Tuple(*var1), Tuple(*var2))
        obj._p = p
        obj._k = k
        obj._var1 = var1
        obj._var2 = var2
        return obj

    def _sympystr(self, printer):
        return printer._print_Function(self)
