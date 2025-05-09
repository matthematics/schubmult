import os
from functools import cache

import symengine.lib.symengine_wrapper as sw
import sympy

# from symengine import SympifyError, sympify
from sympy.core._print_helpers import Printable

os.environ["USE_SYMENGINE"] = "1"
from sympy import sstr
from sympy.core.backend import *
from sympy.core.expr import Expr

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

    @property
    def args(self):
        return self._obj.args

    def __getattr__(self, attr):
        # print(f"{attr}")
        return getattr(self._obj, attr)

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
        return self.obj.args

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
        return self._obj.expand(*args, **kwargs)

    def __repr__(self):
        return self._obj.__repr__()

    def __str__(self):
        return sstr(self._obj)


class SymengineExprClass(type):
    @property
    def __sympyclass__(cls):
        return cls._sympyclass


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

    is_composite: bool | None
    is_noninteger: bool | None
    is_extended_positive: bool | None
    is_negative: bool | None
    is_complex: bool | None
    is_extended_nonpositive: bool | None
    is_integer: bool | None
    is_positive: bool | None
    is_rational: bool | None
    is_extended_nonnegative: bool | None
    is_infinite: bool | None
    is_extended_negative: bool | None
    is_extended_real: bool | None
    is_finite: bool | None
    is_polar: bool | None
    is_imaginary: bool | None
    is_transcendental: bool | None
    is_extended_nonzero: bool | None
    is_nonzero: bool | None
    is_odd: bool | None
    is_algebraic: bool | None
    is_prime: bool | None
    is_commutative: bool | None
    is_nonnegative: bool | None
    is_nonpositive: bool | None
    is_irrational: bool | None
    is_real: bool | None
    is_zero: bool | None
    is_even: bool | None

    # _sympyclass = SympyExpr

    def __new__(cls, *args):
        obj = sw.Symbol.__new__(cls)
        obj._base_args = args
        # print("woo")
        obj._sympyclass = SympyExpr
        return obj

    def __init__(self, *args):
        super().__init__(self, *args, store_pickle=True)
        self._sympy_obj = SympyExpr(self)

    def _sympy_(self):
        return self._sympy_obj

    def __hash__(self):
        return hash(self.args)

    def encode(self, *args):
        from sympy.printing.str import sstr

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
