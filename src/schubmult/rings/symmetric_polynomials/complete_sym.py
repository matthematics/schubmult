from functools import cache

from sympy import Add, Dict, Function, S, sympify

from schubmult.rings.poly_lib import elem_sym_poly
from schubmult.rings.symmetric_polynomials.elem_sym import ElemSym
from schubmult.utils.logging import get_logger

logger = get_logger(__name__)

class h(Function):
    is_commutative = True
    is_Atom = False
    is_polynomial = True
    is_Function = True
    is_nonzero = True

    # def _symengine_(self):
    #     return SymPolyWrap(self, self.args, self.__class__, sw.PyModule(self.__module__))

    def __new__(cls, p, k, *args):
        if hasattr(args[0], "__iter__"):
            return CompleteSym.__xnew_cached__(cls, p, k, tuple(args[0]), tuple(args[1]))
        return CompleteSym.__xnew_cached__(cls, p, k, tuple(args[:k]), tuple(args[k:2*k + p - 1]))

    @staticmethod
    @cache
    def __xnew_cached__(_class, p, k, var1, var2):
        return CompleteSym.__xnew__(_class, p, k, var1, var2)

    def to_elem_sym(self):
        return S.NegativeOne * (self._p % 2) * self._under_elem

    @staticmethod
    def __xnew__(_class, p, k, var1, var2):
        cont_obj = ElemSym(p, k + p - 1, var2, var1)
        if not isinstance(cont_obj, ElemSym):
            return cont_obj
        obj = Function.__new__(_class, cont_obj.args[0], cont_obj.args[1] + 1 - cont_obj.args[0], *cont_obj._coeffvars, *cont_obj._genvars)
        obj._under_elem = cont_obj
        obj._p = obj.args[0]
        obj._k = obj.args[1]
        obj._genvars = cont_obj._coeffvars
        obj._coeffvars = cont_obj._genvars
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
        return self._genvars

    @property
    def coeffvars(self):
        return self._coeffvars

    def _eval_expand_func(self, *args, **kwargs):  # noqa: ARG002
        return sympify(elem_sym_poly(self._under_elem._p, self._under_elem._k, [-x for x in self._under_elem.genvars], [-y for y in self._under_elem.coeffvars]))

    # @property
    # def func(self):
    #     def h(*args):
    #         return self.__class__(*args)
    #     return h

    def _eval_subs(self, *rule):
        rule = Dict(rule)
        new_args = [*self.args]
        for i, arg in enumerate(self.args[2:]):
            new_args[i+2] = arg.xreplace(rule)
        return self.func(*new_args)

    def xreplace(self, rule):
        rule = Dict(rule)
        new_args = [*self.args]
        for i, arg in enumerate(self.args[2:]):
            new_args[i+2] = arg.xreplace(rule)
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

CompleteSym = h
