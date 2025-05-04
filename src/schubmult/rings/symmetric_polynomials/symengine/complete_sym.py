from functools import cache

import sympy
from symengine import Add
from symengine.lib.symengine_wrapper import S, Symbol, sympify
from sympy import Dict, Function, Integer, StrPrinter, Tuple
from sympy.printing.defaults import DefaultPrinting

import schubmult.rings.symmetric_polynomials.sympy.complete_sym as symp
from schubmult.rings.poly_lib import elem_sym_poly
from schubmult.utils.logging import get_logger

from .elem_sym import ElemSym

logger = get_logger(__name__)


class CompleteSym(Symbol, DefaultPrinting):
    is_commutative = True
    is_Atom = False
    is_polynomial = True
    is_Function = True
    is_nonzero = True

    def __new__(cls, p, k, var1, var2):
        return CompleteSym.__xnew_cached__(cls, p, k, tuple(var1), tuple(var2))

    def has_free(self, *x):
        if any(a in self.free_symbols or sympify(a) in self.free_symbols for a in x):
            return True
        return False
    
    @property
    def args(self):
        return(
            Integer(self._p),
            Integer(self._k),
            Tuple(*sorted(self._genvars, key=lambda x: sympy.sympify(x).sort_key())),
            Tuple(*sorted(self._coeffvars, key=lambda x: sympy.sympify(x).sort_key())))
    
    @staticmethod
    @cache
    def __xnew_cached__(_class, p, k, var1, var2):
        return CompleteSym.__xnew__(_class, p, k, var1, var2)

    def to_elem_sym(self):
        return S.NegativeOne * (self._p % 2) * self._under_elem
    
    def __init__(self, *args, **kwargs):
        Symbol.__init__(self, sympy.sstr(self))

    @staticmethod
    def __xnew__(_class, p, k, var1, var2):
        cont_obj = ElemSym(p, k + p - 1, var2, var1)
        if not isinstance(cont_obj, ElemSym):
            return cont_obj
        var1 = var1[:k]
        var2 = var2[:p+k-1]
        name = StrPrinter()._print_Function(Function("h")(p,k,Tuple(*sorted(var1, key=lambda x: sympy.sympify(x).sort_key())),Tuple(*sorted(var2, key=lambda x: sympy.sympify(x).sort_key()))))
        obj = Symbol.__new__(
            _class,
            name,
        )
        obj._p = p
        obj._k = k
        obj._genvars = var1
        obj._coeffvars = var2
        return obj

    @classmethod
    def from_elem_sym(cls, elem):
        return cls(elem._p, elem._k + 1 - elem._p, elem.coeffvars, elem.genvars)

    @property
    def free_symbols(self):
        return set(self.genvars).union(set(self.coeffvars))

    def split_out_vars(self, vars1, vars2=None):  # noqa: ARG002
        if not all(v in self.genvars for v in vars1):
            return self
        first_vars = [sympify(v) for v in vars1 if v in self.genvars]
        second_vars = [a for a in self.genvars if a not in vars1]
        k1 = len(first_vars)
        k2 = len(second_vars)
        if k1 + k2 != self._k:
            raise Exception

        return Add(*[CompleteSym(i, k1, first_vars, self.coeffvars[: k1 + i - 1]) * CompleteSym(self._p - i, k2, second_vars, self.coeffvars[k1 + i :]) for i in range(self._p + 1)])

    @property
    def degree(self):
        return self._p

    @property
    def numvars(self):
        return self._k

    @property
    def genvars(self):
        return tuple(self.args[2])

    @property
    def coeffvars(self):
        return tuple(self.args[3])

    def _eval_expand_func(self, *args, **kwargs):  # noqa: ARG002
        return sympify(elem_sym_poly(self._under_elem._p, self._under_elem._k, [-x for x in self._under_elem.genvars], [-y for y in self._under_elem.coeffvars]))

    @property
    def func(self):
        def h(*args):
            return self.__class__(*args)

        return h

    def _eval_subs(self, *rule):
        rule = Dict(rule)
        new_args = [*self.args]
        new_args[2] = [*new_args[2]]
        new_args[3] = [*new_args[3]]
        for i, arg in enumerate(self.args[2]):
            new_args[2][i] = arg.subs(rule)
        for i, arg in enumerate(self.args[3]):
            new_args[3][i] = arg.subs(rule)
        return self.func(*new_args)

    def xreplace(self, rule):
        rule = Dict(rule)
        new_args = [*self.args]
        new_args[2] = [*new_args[2]]
        new_args[3] = [*new_args[3]]
        for i, arg in enumerate(self.args[2]):
            new_args[2][i] = arg.xreplace(rule)
        for i, arg in enumerate(self.args[3]):
            new_args[3][i] = arg.xreplace(rule)
        return self.func(*new_args)

    def divide_out_diff(self, v1, v2):
        new_obj = self._under_elem.divide_out_diff(v1, v2)
        return new_obj.replace(ElemSym, lambda x: CompleteSym.from_elem_sym(ElemSym(*x)))

    @staticmethod
    def from_expr_elem_sym(expr):
        return expr.replace(ElemSym, lambda x: CompleteSym.from_elem_sym(ElemSym(*x)))

    def _eval_div_diff(self, v1, v2):
        return CompleteSym.from_expr_elem_sym(self._under_elem.div_diff(v1, v2))

    def _eval_divide_out_diff(self, v1, v2):
        return CompleteSym.from_expr_elem_sym(self._under_elem.divide_out_diff(v1, v2))

    def div_diff(self, v1, v2):
        return CompleteSym.from_expr_elem_sym(self._under_elem.div_diff(v1, v2))

    def _sympystr(self, printer):
        return printer._print_Function(self)
    
    def __str__(self):
        return self._sympy_().__str__() #sympy.sstr(self._obj)

    def _sympy_(self):
        return symp.CompleteSym(*self.args)

    def _sympystr(self, printer):
        return printer._print_Function(self._sympy_())