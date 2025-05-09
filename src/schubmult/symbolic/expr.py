import os
from functools import cache

import symengine.lib.symengine_wrapper as sw
import sympy

# from symengine import SympifyError, sympify
from sympy.core._print_helpers import Printable
from sympy.core.expr import Expr

from .base_printing import sstr

# class SympyExprClass(type):
#     @property
#     def __symengineclass__(cls):
#         return cls._symengineclass

#     def __repr__(cls):
#         return repr(cls._symengineclass)


class SympyExpr(Expr):
    def __new__(cls, _obj):
        obj = Expr.__new__(cls, *_obj.args)
        obj._obj = _obj
        obj.__dict__.update(_obj.__dict__)
        return obj

    def subs(self, *args, **kwargs):
        if hasattr(self._obj, "subs"):
            return self._obj.subs(*args, **kwargs)
        return super().subs(*args, **kwargs)

    def xreplace(self, *args, **kwargs):
        if hasattr(self._obj, "xreplace"):
            return self._obj.xreplace(*args, **kwargs)
        return super().xreplace(*args, **kwargs)

    def _eval_subs(self, *args, **kwargs):
        if hasattr(self._obj, "_eval_subs"):
            return self._obj.subs(*args, **kwargs)
        return super()._eval_subs(*args, **kwargs)

    def __hash__(self):
        return hash(self.args)

    @property
    def is_number(self):
        return False

    def _symengine_(self):
        # print("pageblitz")
        return self._obj

    def _sympystr(self, printer):
        return printer.doprint(self._obj)

    # @property
    # def args(self):
    #     return [sympy.sympify(arg) for arg in self._obj.args]

    def __getattr__(self, attr):
        # print(f"{self=} {attr}")
        pang = getattr(self._obj, attr)
        # print(f"{pang=}")
        return pang

    def _sympy_(self):
        return self

    # need to explicitly forward overriden methods by sympy expr

    def simplify(self, **kwargs) -> Expr:
        if hasattr(self._obj, "simplify"):
            return sympy.sympify(self._obj.simplify(**kwargs))
        return super().simplify(**kwargs)

    @cache
    def sort_key(self, order=None):
        if hasattr(self._obj, "sort_key"):
            return self._obj.sort_key(order)
        return super().sort_key(order)

    def _hashable_content(self):
        return self._obj.args

    def equals(self, other, failing_expression=False):
        if hasattr(self._obj, "equals"):
            return self._obj.equals(other)
        return super().equals(other, failing_expression=failing_expression)

    def as_ordered_factors(self, order=None):
        if hasattr(self._obj, "as_ordered_factors"):
            return self._obj.as_ordered_factors(order)
        return super().as_ordered_factors(order)

    # def as_poly(self, *gens, **args):

    def as_ordered_terms(self, order=None):
        if hasattr(self._obj, "as_ordered_terms"):
            return self._obj.as_ordered_terms(order)
        return super().as_ordered_terms(order)

    def as_terms(self):
        if hasattr(self._obj, "as_terms"):
            return self._obj.as_terms()
        return super().as_terms()

    def as_powers_dict(self):
        if hasattr(self._obj, "as_powers_dict"):
            return self._obj.as_powers_dict()
        return super().as_powers_dict()

    def as_coefficients_dict(self, *syms):
        if hasattr(self._obj, "as_coefficients_dict"):
            return self._obj.as_coefficients_dict(*syms)
        return super().as_coefficients_dict(*syms)

    def as_base_exp(self):
        if hasattr(self._obj, "as_base_exp"):
            return self._obj.as_base_exp()
        return super().as_base_exp()

    def func(self, *args):
        return sympy.sympify(self._obj.func(*args))

    # def could_extract_minus_sign(self) -> bool:
    #     """Return True if self has -1 as a leading factor or has
    #     more literal negative signs than positive signs in a sum,
    #     otherwise False.

    #     Examples
    #     ========

    #     >>> from sympy.abc import x, y
    #     >>> e = x - y
    #     >>> {i.could_extract_minus_sign() for i in (e, -e)}
    #     {False, True}

    #     Though the ``y - x`` is considered like ``-(x - y)``, since it
    #     is in a product without a leading factor of -1, the result is
    #     false below:

    #     >>> (x*(y - x)).could_extract_minus_sign()
    #     False

    #     To put something in canonical form wrt to sign, use `signsimp`:

    #     >>> from sympy import signsimp
    #     >>> signsimp(x*(y - x))
    #     -x*(x - y)
    #     >>> _.could_extract_minus_sign()
    #     True
    #     """
    #     return False

    # @staticmethod
    # def _expand_hint(expr, hint, deep=True, **hints):
    #     """
    #     Helper for ``expand()``.  Recursively calls ``expr._eval_expand_hint()``.

    #     Returns ``(expr, hit)``, where expr is the (possibly) expanded
    #     ``expr`` and ``hit`` is ``True`` if ``expr`` was truly expanded and
    #     ``False`` otherwise.
    #     """
    #     hit = False
    #     # XXX: Hack to support non-Basic args
    #     #              |
    #     #              V
    #     if deep and getattr(expr, 'args', ()) and not expr.is_Atom:
    #         sargs = []
    #         for arg in expr.args:
    #             arg, arghit = Expr._expand_hint(arg, hint, **hints)
    #             hit |= arghit
    #             sargs.append(arg)

    #         if hit:
    #             expr = expr.func(*sargs)

    #     if hasattr(expr, hint):
    #         newexpr = getattr(expr, hint)(**hints)
    #         if newexpr != expr:
    #             return (newexpr, True)

    #     return (expr, hit)

    def expand(self, *args, **kwargs):
        try:
            return self._obj.expand(*args, **kwargs)
        except Exception:
            return super().expand(*args, **kwargs)

    def __repr__(self):
        return self._obj.__repr__()

    def __str__(self):
        return sstr(self._obj)

    def compare(self, other):
        # print("bongfunket")
        if hasattr(self._obj, "compare"):
            return self._obj.compare(other)
        return super().compare(other)
    
    def has(self, *args):
        if hasattr(self._obj, "has"):
            return self._obj.has(*args)
        return super().has(*args)

    def has_symbol(self, sym):
        if hasattr(self._obj, "has_symbol"):
            return self._obj.has_symbol(sym)
        return super().has_symbol(sym)


class SymengineExprClass(type):
    @property
    def __class__(self):
        return sw.Expr


class SymengineExpr(sw.Symbol, Printable, metaclass=SymengineExprClass):
    _op_priority = 800000

    is_number = False
    is_Atom = False
    is_Symbol = False
    is_symbol = False
    is_Indexed = False
    is_Dummy = False
    is_Wild = False
    is_Function = False
    is_Add = False
    is_Mul = False
    is_Pow = False
    is_Number = False
    is_Float = False
    is_Rational = False
    is_Integer = False
    is_NumberSymbol = False
    is_Order = False
    is_Derivative = False
    is_Piecewise = False
    is_Poly = False
    is_AlgebraicNumber = False
    is_Relational = False
    is_Equality = False
    is_Boolean = False
    is_Not = False
    is_Matrix = False
    is_Vector = False
    is_Point = False
    is_MatAdd = False
    is_MatMul = False

    # _sympyclass = SympyExpr

    def __str__(self):
        return sstr(self)

    def __getattr__(self, attr):
        # print(f"{self=} {attr=}")
        raise AttributeError

    def __new__(cls, *args):
        obj = sw.Symbol.__new__(cls)
        # print("{args=}")
        obj._base_args = args
        # print("woo")
        obj._sympyclass = SympyExpr
        return obj

    def __init__(self, *args):
        # print([type(arg) for arg in args])
        # print(f"{type(self)=}")
        sw.Symbol.__init__(self, self, *args, store_pickle=False)
        # self._sympy_obj = SympyExpr(self)

    def _sympy_(self):
        return SympyExpr(self)

    def __hash__(self):
        return hash(self.args)

    def encode(self, *args):
        return sstr(self).encode(*args)

    @property
    def args(self):
        return self._base_args

    def __setstate__(self, state):
        self.__dict__.update(state)

    def __reduce__(self):
        return (self.__class__, self.args, self.__dict__)

    def _sympystr(self, printer):
        return printer._print(self.args)

    def _sympyrepr(self, printer):
        return printer._print(self.args)

    @property
    def func(self):
        return self.__class__
