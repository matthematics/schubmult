from functools import cache

import symengine.lib.symengine_wrapper as sw
import sympy
from symengine import SympifyError, sympify
from sympy.core._print_helpers import Printable

from schubmult.rings.variables import CustomGeneratingSet, GeneratingSet, ZeroGeneratingSet
from schubmult.utils.perm_utils import add_perm_dict

NoneVar = 1e10
ZeroVar = 0


class NotEnoughGeneratorsError(ValueError):
    pass


@cache
def poly_ring(v: str):
    if v == ZeroVar:
        return ZeroGeneratingSet(tuple([sympify(0) for i in range(100)]))
    if v == NoneVar:
        return ZeroGeneratingSet(tuple([sympify(0) for i in range(100)]))
    return GeneratingSet(str(v))


def _mul_schub_dicts(dict1, dict2, basis1, basis2, best_effort_positive=True):
    this_dict = {}
    for k, v in dict2.items():
        for kd, vd in dict1.items():
            did_positive = False
            to_mul = v * vd
            if best_effort_positive:
                try:
                    this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis1.cached_positive_product(kd, k, basis2).items()})
                    did_positive = True
                except Exception:
                    did_positive = False
            if not did_positive:
                this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis1.cached_product(kd, k, basis2).items()})
    return this_dict


def _tensor_product_of_dicts(d1, d2):
    ret_dict = {}
    for k1, v1 in d1.items():
        this_dict = {}
        try:
            v1 = sympify(v1)
        except SympifyError:
            v1 = sympy.sympify(v1)
        for k2, v2 in d2.items():
            try:
                v2 = sympify(v2)
            except SympifyError:
                v2 = sympy.sympify(v2)
                v1 = sympy.sympify(v1)
            if isinstance(k1, tuple):
                this_dict[(*k1, k2)] = v1 * v2
            else:
                this_dict[(k1, k2)] = v1 * v2
        ret_dict = add_perm_dict(ret_dict, this_dict)
    return ret_dict


import os

os.environ["USE_SYMENGINE"] = "1"
from sympy import sstr
from sympy.core.backend import *
from sympy.core.expr import Expr


class SympyExprClass(type):
    @property
    def __symengineclass__(cls):
        return cls._symengineclass

    def __repr__(cls):
        return repr(cls._symengineclass)


class SympyExpr(Expr, metaclass=SympyExprClass):

    def __new__(cls, _obj):
        obj = Expr.__new__(cls)
        obj._obj = _obj
        obj.__dict__.update(_obj.__dict__)
        return obj

    def __hash__(self):
        return hash(self.args)

    @property
    def is_number(self):
        return False
    
    # def __getattr__(self, attr):
    #     if attr == "_symengine_":
    #         return self._symengine_
    #     return getattr(self._obj, attr)

    def _symengine_(self):
        return self._obj

    def _sympystr(self, printer):
        return printer.doprint(self._obj)

    @property
    def args(self):
        return self._obj.args

    @property
    def func(self):
        def func(*bob):
            return self.__class__(self._obj.func(*bob))
        return func

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

    #_sympyclass = SympyExpr
    js_sympy = False
    def __new__(cls, *args):
        obj = sw.Symbol.__new__(cls)
        obj._base_args = args
        class Me(SympyExpr):
            _symengineclass = cls
            @property
            def __name__(self):
                return Me._symengineclass.__name__

        obj._sympyclass = Me
        obj._obj = Me.__new__(obj._sympyclass, obj)
        return obj

    def __init__(self, *args):
        super().__init__(self, *args, store_pickle=True)

    def _sympy_(self):
        return self._obj

    def __hash__(self):
        return hash(self.args)

    def encode(self, *args):
        from sympy.printing.str import sstr

        return sstr(self).encode(*args)

    @property
    def args(self):
        return self._base_args

    def __reduce__(self):
        return (self.__class__, self.args)

    def _sympystr(self, printer):
        return printer._print(self.args)

    def _sympyrepr(self, printer):
        return printer._print(self.args)

    @property
    def func(self):
        return self.__class__
