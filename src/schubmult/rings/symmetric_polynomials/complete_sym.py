from functools import cache

from schubmult.rings.poly_lib import elem_sym_poly
from schubmult.rings.symmetric_polynomials.elem_sym import FactorialElemSym
from schubmult.rings.variables import NotEnoughGeneratorsError, ZeroGeneratingSet
from schubmult.symbolic import Function, Integer, S, sympify, sympify_sympy, sympy_Add
from schubmult.utils.logging import get_logger

logger = get_logger(__name__)


class CompleteSym_base(Function):
    def _eval_subs(self, rule):
        new_args = [*self.args]
        for i, arg in enumerate(self.args[2:]):
            new_args[i + 2] = sympify(arg).xreplace(rule)
        return self.func(*new_args)

    def xreplace(self, rule):
        new_args = [*self.args]
        for i, arg in enumerate(self.args[2:]):
            new_args[i + 2] = sympify(arg).xreplace(rule)
        return self.func(*new_args)

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

    def _sympystr(self, printer):
        return printer._print_Function(self)


class H(CompleteSym_base):
    is_commutative = True
    is_Atom = False
    is_polynomial = True
    is_Function = True
    is_nonzero = True

    def __new__(cls, p, k, *args):
        p = int(p)
        k = int(k)
        if hasattr(args[0], "__iter__"):
            return FactorialCompleteSym.__xnew_cached__(cls, int(p), int(k), tuple(args[0]), tuple(args[1]))
        return FactorialCompleteSym.__xnew_cached__(cls, int(p), int(k), tuple(args[: int(k)]), tuple(args[k:]))

    @staticmethod
    @cache
    def __xnew_cached__(_class, p, k, var1, var2):
        return FactorialCompleteSym.__xnew__(_class, p, k, var1, var2)

    def to_elem_sym(self):
        return S.NegativeOne * (self._p % 2) * self._under_elem

    @staticmethod
    def __xnew__(_class, p, k, var1, var2):
        cont_obj = FactorialElemSym(p, k + p - 1, var2, var1)
        if not isinstance(cont_obj, FactorialElemSym):
            return cont_obj
        obj = CompleteSym_base.__new__(_class, cont_obj.args[0], cont_obj.args[1] + 1 - cont_obj.args[0], *cont_obj._coeffvars, *cont_obj._genvars)
        obj._under_elem = cont_obj
        obj._p = int(obj.args[0])
        obj._k = int(obj.args[1])
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

        return sympy_Add(
            *[FactorialCompleteSym(i, k1, first_vars, self.coeffvars[: k1 + i - 1]) * FactorialCompleteSym(self._p - i, k2, second_vars, self.coeffvars[k1 + i :]) for i in range(self._p + 1)]
        )

    def _eval_expand_func(self, *args, **kwargs):  # noqa: ARG002
        return sympify_sympy(elem_sym_poly(self._under_elem._p, self._under_elem._k, [-x for x in self._under_elem.genvars], [-y for y in self._under_elem.coeffvars]))

    @property
    def func(self):
        def H(*args):
            return self.__class__(*args)

        return H

    def divide_out_diff(self, v1, v2):
        new_obj = self._under_elem.divide_out_diff(v1, v2)
        return new_obj.replace(FactorialElemSym, lambda x: FactorialCompleteSym.from_elem_sym(FactorialElemSym(*x)))

    @staticmethod
    def from_expr_elem_sym(expr):
        return expr.replace(FactorialElemSym, lambda x: FactorialCompleteSym.from_elem_sym(FactorialElemSym(*x)))

    def _eval_div_diff(self, v1, v2):
        return FactorialCompleteSym.from_expr_elem_sym(self._under_elem.div_diff(v1, v2))

    def _eval_divide_out_diff(self, v1, v2):
        return FactorialCompleteSym.from_expr_elem_sym(self._under_elem.divide_out_diff(v1, v2))

    def div_diff(self, v1, v2):
        return FactorialCompleteSym.from_expr_elem_sym(self._under_elem.div_diff(v1, v2))


class h(CompleteSym_base):
    is_commutative = True
    is_Atom = False
    is_polynomial = True
    is_Function = True
    is_nonzero = True

    # _sympyclass = _h
    # def _symengine_(self):
    #     return SymPolyWrap(self, self.args, self.__class__, sw.PyModule(self.__module__))

    def __new__(cls, p, k, *args):
        p = int(p)
        k = int(k)
        if hasattr(args[0], "__iter__"):
            return CompleteSym.__xnew_cached__(cls, int(p), int(k), tuple(args[0]))
        return CompleteSym.__xnew_cached__(cls, int(p), int(k), tuple(args))

    # def __hash__(self):
    #     return hash(self.args, "bonbon")

    @staticmethod
    @cache
    def __xnew_cached__(_class, p, k, var1):
        return CompleteSym.__xnew__(_class, p, k, var1)

    # def to_elem_sym(self):
    #     return S.NegativeOne * (self._p % 2) * self._under_elem

    @staticmethod
    def __xnew__(_class, p, k, var1):
        if len(var1) < k:
            raise NotEnoughGeneratorsError
        var1 = var1[:k]
        obj = CompleteSym_base.__new__(_class, Integer(p), Integer(k), *var1)
        obj._p = int(obj.args[0])
        obj._k = int(obj.args[1])
        obj._genvars = tuple(var1)
        return obj

    @property
    def free_symbols(self):
        return set(self.genvars)

    @property
    def coeffvars(self):
        return ZeroGeneratingSet()

    def _eval_expand_func(self, *args, **kwargs):  # noqa: ARG002
        return sympify(elem_sym_poly(self._p, self._p + self._k - 1, self.coeffvars, [-y for y in self.genvars]))

    @property
    def func(self):
        def h(*args):
            return self.__class__(*args)

        return h


FactorialCompleteSym = H
CompleteSym = h
