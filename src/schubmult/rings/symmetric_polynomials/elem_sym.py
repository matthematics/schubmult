from functools import cache

import symengine.lib.symengine_wrapper as sw
from sympy import Dict, FiniteSet, Integer, S, Symbol, Tuple, sympify
from sympy.core.function import Function

# import schubmult.rings.symmetric_polynomials.symengine.elem_sym as syme
from schubmult.rings.poly_lib import elem_sym_poly
from schubmult.utils.logging import get_logger
from schubmult.utils.ring_utils import NotEnoughGeneratorsError

logger = get_logger(__name__)


class E(Function):
    is_commutative = True
    is_Atom = False
    is_polynomial = True
    is_Function = True
    is_nonzero = True

    def __new__(cls, p, k, *args):
        if hasattr(args[0], "__iter__"):
            return ElemSym.__xnew_cached__(cls, p, k, tuple(args[0]), tuple(args[1]))
        return ElemSym.__xnew_cached__(cls, p, k, tuple(args[:k]), tuple(args[k:2*k + 1 - p]))

    @staticmethod
    @cache
    def __xnew_cached__(_class, p, k, var1, var2):
        return ElemSym.__xnew__(_class, p, k, var1, var2)

    @staticmethod
    def __xnew__(_class, p, k, var1, var2):
        if p > k or k < 0:
            return S.Zero
        if p == 0:
            return S.One
        var1 = var1[:k]
        var2 = var2[: k + 1 - p]
        for i, v in enumerate(var1):
            if v in var2:
                j = var2.index(v)
                return ElemSym.__new__(_class, p, k - 1, [*var1[:i], *var1[i + 1 :]], [*var2[:j], *var2[j + 1 :]])
        if len(var1) < k:
            raise NotEnoughGeneratorsError(f"{k} passed as number of variables but only {len(var1)} given")
        if len(var2) < k + 1 - p:
            raise NotEnoughGeneratorsError(f"{k} passed as number of variables and degree is {p} but only {len(var2)} coefficient variables given. {k + 1 - p} coefficient variables are needed.")
        var1 = tuple(sorted(var1,key=lambda x: sympify(x).sort_key()))
        var2 = tuple(sorted(var2,key=lambda x: sympify(x).sort_key()))
        obj = Function.__new__(
            _class,
            Integer(p),
            Integer(k),
            *var1,
            *var2,
        )
        # if len(obj.args[2]) < k:
        #     raise ValueError("Duplicate genvar arguments")
        # if len(obj.args[3]) < k + 1 - p:
        #     raise ValueError("Duplicate coeffvar arguments")
        obj._p = p
        obj._k = k
        obj._genvars = var1
        obj._coeffvars = var2

        return obj


    def _symengine_(self):
        return sw.PyFunction(self, (self.args[0],self.args[1],*self.genvars,*self.coeffvars), self.__class__, sw.PyModule(self.__module__))

    @property
    def free_symbols(self):
        return set(self.genvars).union(set(self.coeffvars))

    def split_out_vars(self, vars1, vars2=None):
        # order of vars2 matters!
        vars1 = [sympify(v) for v in vars1]
        if vars2 is None:
            vars2 = [*self.coeffvars]
        else:
            vars2 = [sympify(v) for v in vars2]
        if not all(v in self.genvars for v in vars1):
            raise NotEnoughGeneratorsError(f"Not all variables {vars1} are in the generating set {self.genvars}")
        newvars2 = [*vars2]
        for v in vars2:
            if v not in self.coeffvars:
                newvars2.remove(v)

        # vars2 = [*newvars2]
        # for v in vars1:
        #     if v not in self.genvars:
        #         new_vars1.remove(v)
        # new_vars2 = [*vars2]
        # for v in vars2:
        #     if v not in self.coeffvars:
        #         new_vars2.remove(v)
        # vars1 = new_vars1
        # vars2 = new_vars2
        ret = S.Zero
        k1 = len(vars1)
        k2 = self._k - k1
        new_genvars1 = [*self.genvars]
        for v in vars1:
            new_genvars1.remove(v)
        # print(len(new_genvars1))
        new_coeffvars2 = [*newvars2]
        new_coeffvars1 = [*newvars2]
        new_coeffvars2.reverse()
        # we have k + 1 - p
        for p1 in range(min(k1 + 1, self._p + 1)):
            p2 = self._p - p1
            # print(f"{p1=}, {p2=}, {k1=}, {k2=}")
            # print(f"{vars1=}, {vars2=}")
            # print(f"{new_coeffvars1=} {new_coeffvars2}")
            try:
                ret += self.func(p1, k1, vars1, new_coeffvars1) * self.func(p2, k2, new_genvars1, new_coeffvars2)
            except NotEnoughGeneratorsError:
                pass
        return ret

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

    @cache
    def _eval_expand_func(self, *args, **kwargs):  # noqa: ARG002
        return sympify(elem_sym_poly(self._p, self._k, self.genvars, self.coeffvars))

    # @property
    # def func(self):
    #     def e(*args):
    #         return self.__class__(*args)
    #     return e

    #     return e

    #     return e
    def _eval_subs(self, *rule):
        # print(f"_eval_subs")
        # print(f"{rule=}")
        # print(f"{self=}")
        rule = Dict(rule)
        new_args = [*self.args]
        new_args = [*self.args]
        for i, arg in enumerate(self.args[2:]):
            new_args[i + 2] = arg.xreplace(rule)
        return self.func(*new_args)

    def xreplace(self, rule):
        # print(f"xreplace")
        # print(f"{rule=}")
        # print(f"{self=}")
        rule = Dict(rule)
        new_args = [*self.args]
        for i, arg in enumerate(self.args[2:]):
            new_args[i + 2] = arg.xreplace(rule)
        return self.func(*new_args)

    def divide_out_diff(self, v1, v2):
        if v1 == v2:
            return S.Zero
        if v1 in self.genvars:
            new_genvars = [*self.genvars]
            new_genvars.remove(v1)
            return self.func(self._p - 1, self._k - 1, new_genvars, self.coeffvars)
        if v1 in self.coeffvars:
            if v2 in self.coeffvars:
                return S.Zero
            if v2 in self.genvars:
                new_genvars = [*self.genvars]
                new_genvars.remove(v2)
                return -self.func(self._p, self._k - 1, new_genvars, self.coeffvars)
            return -self.func(self._p - 1, self._k, self.genvars, [*self.coeffvars, v2])
        return S.Zero

    def _eval_div_diff(self, v1, v2):
        return self.div_diff(v1, v2)

    def _eval_divide_out_diff(self, v1, v2):
        return self.divide_out_diff(v1, v2)

    def div_diff(self, v1, v2):
        if v1 == v2:
            return S.Zero
        if v1 in self.genvars:
            new_genvars = [*self.genvars]
            new_genvars.remove(v1)
            return self.func(self._p - 1, self._k - 1, new_genvars, self.coeffvars)
        if v1 in self.coeffvars:
            if v2 in self.coeffvars:
                return S.Zero
            if v2 in self.genvars:
                new_genvars = [*self.genvars]
                new_genvars.remove(v2)
                return -self.func(self._p, self._k - 1, new_genvars, self.coeffvars)
            return -self.func(self._p - 1, self._k, self.genvars, [*self.coeffvars, v2])
        if v2 in self.genvars:
            new_genvars = [*self.genvars]
            new_genvars.remove(v2)
            return -self.func(self._p - 1, self._k - 1, new_genvars, self.coeffvars)
        if v2 in self.coeffvars:
            if v1 in self.coeffvars:
                return S.Zero
            if v1 in self.genvars:
                new_genvars = [*self.genvars]
                new_genvars.remove(v1)
                return self.func(self._p, self._k - 1, new_genvars, self.coeffvars)
            return self.func(self._p - 1, self._k, self.genvars, [*self.coeffvars, v1])
        return S.Zero

    def sign_of_pair(self, v1, v2):
        if v1 == v2:
            return S.Zero
        if v1 in self.genvars:
            return S.One
        if v1 in self.coeffvars:
            if v2 in self.coeffvars:
                return S.Zero
            if v2 in self.genvars:
                return None
            return S.NegativeOne
        return S.Zero

    def _sympystr(self, printer):
        # return printer._print_Function(self)
        return printer._print_Function(self)

    def pull_out_vars(self, var1, var2, min_degree=1):
        if self._p < min_degree:
            return self
        if var1 in self.genvars and var2 in self.coeffvars:
            return self.xreplace({var1: var2}) + ElemSym(1, 1, [var1], [var2]) * self.divide_out_diff(var1, var2)
        return self


class e(Function):
    is_commutative = True
    is_Atom = False
    is_polynomial = True
    is_Function = True
    is_nonzero = True

    def __new__(cls, p, k, *args):
        if hasattr(args[0], "__iter__"):
            return ElemSym.__xnew_cached__(cls, p, k, tuple(args[0]))
        return ElemSym.__xnew_cached__(cls, p, k, tuple(args))

    @staticmethod
    @cache
    def __xnew_cached__(_class, p, k, var1):
        return ElemSym.__xnew__(_class, p, k, var1)

    @staticmethod
    def __xnew__(_class, p, k, var1, var2):
        if p > k or k < 0:
            return S.Zero
        if p == 0:
            return S.One
        var1 = var1[:k]
        if len(var1) < k:
            raise NotEnoughGeneratorsError(f"{k} passed as number of variables but only {len(var1)} given")
        var1 = tuple(sorted(var1,key=lambda x: sympify(x).sort_key()))
        obj = Function.__new__(
            _class,
            Integer(p),
            Integer(k),
            *var1,
        )
        # if len(obj.args[2]) < k:
        #     raise ValueError("Duplicate genvar arguments")
        # if len(obj.args[3]) < k + 1 - p:
        #     raise ValueError("Duplicate coeffvar arguments")
        obj._p = p
        obj._k = k
        obj._genvars = var1

        return obj


    def _symengine_(self):
        return sw.PyFunction(self, (self.args[0],self.args[1],*self.genvars), self.__class__, sw.PyModule(self.__module__))

    @property
    def free_symbols(self):
        return set(self.genvars).union(set(self.coeffvars))

    def split_out_vars(self, vars1, vars2=None):
        # order of vars2 matters!
        vars1 = [sympify(v) for v in vars1]
        if vars2 is not None:
            raise ValueError(f"There are no coeffvars for {self}")

        if not all(v in self.genvars for v in vars1):
            raise NotEnoughGeneratorsError(f"Not all variables {vars1} are in the generating set {self.genvars}")

        ret = S.Zero
        k1 = len(vars1)
        k2 = self._k - k1
        new_genvars1 = [*self.genvars]
        for v in vars1:
            new_genvars1.remove(v)
        # print(len(new_genvars1))
        for p1 in range(min(k1 + 1, self._p + 1)):
            p2 = self._p - p1
            # print(f"{p1=}, {p2=}, {k1=}, {k2=}")
            # print(f"{vars1=}, {vars2=}")
            # print(f"{new_coeffvars1=} {new_coeffvars2}")
            try:
                ret += self.func(p1, k1, vars1) * self.func(p2, k2, new_genvars1)
            except NotEnoughGeneratorsError:
                pass
        return ret

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
        return None

    @cache
    def _eval_expand_func(self, *args, **kwargs):  # noqa: ARG002
        return sympify(elem_sym_poly(self._p, self._k, self.genvars, [0 for i in range(30)]))


    def _eval_subs(self, *rule):
        # print(f"_eval_subs")
        # print(f"{rule=}")
        # print(f"{self=}")
        rule = Dict(rule)
        new_args = [*self.args]
        new_args = [*self.args]
        for i, arg in enumerate(self.args[2:]):
            new_args[i + 2] = arg.xreplace(rule)
        return self.func(*new_args)

    def xreplace(self, rule):
        # print(f"xreplace")
        # print(f"{rule=}")
        # print(f"{self=}")
        rule = Dict(rule)
        new_args = [*self.args]
        for i, arg in enumerate(self.args[2:]):
            new_args[i + 2] = arg.xreplace(rule)
        return self.func(*new_args)

    # let vars be a set
    def divide_out_diff(self, v1, v2):
        if v1 == v2:
            return S.Zero
        if v1 in self.genvars:
            new_genvars = [*self.genvars]
            new_genvars.remove(v1)
            return self.func(self._p - 1, self._k - 1, new_genvars)
        return S.Zero

    def _eval_div_diff(self, v1, v2):
        return self.div_diff(v1, v2)

    def _eval_divide_out_diff(self, v1, v2):
        return self.divide_out_diff(v1, v2)

    def div_diff(self, v1, v2):
        if v1 == v2:
            return S.Zero
        if v1 in self.genvars:
            new_genvars = [*self.genvars]
            new_genvars.remove(v1)
            return self.func(self._p - 1, self._k - 1, new_genvars)

        if v2 in self.genvars:
            new_genvars = [*self.genvars]
            new_genvars.remove(v2)
            return -self.func(self._p - 1, self._k - 1, new_genvars)
        return S.Zero

    def sign_of_pair(self, v1, v2):
        if v1 == v2:
            return S.Zero
        if v1 in self.genvars:
            return S.One
        return S.Zero

    def _sympystr(self, printer):
        # return printer._print_Function(self)
        return printer._print_Function(self)

    def pull_out_vars(self, var1, var2, min_degree=1):
        if self._p < min_degree:
            return self
        if var1 in self.genvars and var2 in self.coeffvars:
            return self.xreplace({var1: var2}) + ElemSym(1, 1, [var1], [var2]) * self.divide_out_diff(var1, var2)
        return self

FactorialElemSym = E

ElemSym = e




