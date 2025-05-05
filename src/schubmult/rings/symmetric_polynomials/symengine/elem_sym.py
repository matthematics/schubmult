import os

os.environ['USE_SYMENGINE'] = "1"



from functools import cache

from symengine.lib.symengine_wrapper import S, Symbol, sympify
from sympy import Dict, Function, Integer, StrPrinter, Tuple, sstr
from sympy.printing.defaults import DefaultPrinting

from schubmult.rings.poly_lib import elem_sym_poly
from schubmult.utils.logging import get_logger
from schubmult.utils.ring_utils import NotEnoughGeneratorsError

logger = get_logger(__name__)


from sympy.core.backend import *


class ElemSym(Symbol, DefaultPrinting):
    is_commutative = True
    is_Atom = True
    is_polynomial = True
    is_Function = True
    is_nonzero = True

    def __new__(cls, p, k, var1, var2):
        return ElemSym.__xnew_cached__(cls, p, k, tuple(var1), tuple(var2))

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
        name = StrPrinter()._print_Function(Function("e")(p,k,*var1, *var2))
        obj = Symbol.__new__(
            _class,
            name,
        )
        obj._p = p
        obj._k = k
        obj._genvars = var1
        obj._coeffvars = var2
        return obj

    def __hash__(self):
        return hash(self.args)

    def __init__(self, *args, **kwargs):
        Symbol.__init__(self, sstr(self))

    @property
    def args(self):
        return(
            Integer(self._p),
            Integer(self._k),
            Tuple(*self._genvars),
            Tuple(*self._coeffvars))

    #def __rmul__(self, other):

    @property
    def free_symbols(self):
        return set(self.genvars).union(set(self.coeffvars))

    def has_free(self, *x):
        if any(a in self.free_symbols or sympify(a) in self.free_symbols for a in x):
            return True
        return False

    def split_out_vars(self, vars1, vars2=None):
        vars1 = [sympify(v) for v in vars1 if v in self.genvars]
        if vars2 is None:
            vars2 = [*self.coeffvars]
        else:
            vars2 = [sympify(v) for v in vars2]
        newvars1 = [*vars1]
        for v in vars1:
            if v not in self.genvars:
                newvars1.remove(v)
        if not len(newvars1):
            return self
        newvars2 = [*vars2]
        for v in vars2:
            if v not in self.coeffvars:
                newvars2.remove(v)

        ret = S.Zero
        k1 = len(newvars1)
        k2 = self._k - k1
        new_genvars1 = [*self.genvars]
        for v in newvars1:
            new_genvars1.remove(v)
        # print(len(new_genvars1))
        new_coeffvars2 = [*newvars2]
        new_coeffvars1 = [*newvars2]
        new_coeffvars2.reverse()
        # we have k + 1 - p
        for p1 in range(min(k1 + 1, self._p + 1)):
            p2 = self._p - p1
            try:
                ret += self.func(p1, k1, newvars1, new_coeffvars1) * self.func(p2, k2, new_genvars1, new_coeffvars2)
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

    @property
    def func(self):
        def e(*args):
            return self.__class__(*args)
        return e

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

    def pull_out_vars(self, var1, var2, min_degree=1):
        if self._p < min_degree:
            return self
        if var1 in self.genvars and var2 in self.coeffvars:
            return self.xreplace({var1: var2}) + ElemSym(1, 1, [var1], [var2]) * self.divide_out_diff(var1, var2)
        return self

    def __str__(self):
        return sstr(self)

    # def _sympy_(self):
    #     return symp.ElemSym(*self.args)

    def _sympystr(self, printer):
        return printer._print_Function(self)


