from functools import cache

from sympy import Integer, S, Tuple, sympify
from sympy.core.expr import Expr

from schubmult.poly_lib.poly_lib import elem_sym_poly

from ._utils import NotEnoughGeneratorsError


class ElemSym(Expr):
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
        obj = Expr.__new__(
            _class,
            Integer(p),
            Integer(k),
            Tuple(*sorted(var1, key=lambda x: int(x.name.split("_")[1]))),
            Tuple(*sorted(var2, key=lambda x: int(x.name.split("_")[1]))),
        )
        obj._p = p
        obj._k = k
        return obj

    def split_out_vars(self, vars1, vars2):
        # order of vars2 matters!
        vars1 = [sympify(v) for v in vars1]
        vars2 = [sympify(v) for v in vars2]
        if not all(v in self.genvars for v in vars1):
            raise NotEnoughGeneratorsError(f"Not all variables {vars1} are in the generating set {self.genvars}")
        if not all(v in self.coeffvars for v in vars2):
            raise NotEnoughGeneratorsError(f"Not all variables {vars2} are in the coefficient set {self.coeffvars}")
        ret = S.Zero
        k2 = len(vars1)
        k1 = self._k - k2
        new_genvars1 = [*self.genvars]
        for v in vars1:
            new_genvars1.remove(v)
        # print(len(new_genvars1))
        new_coeffvars2 = [*vars2]
        new_coeffvars1 = [*vars2]
        new_coeffvars2.reverse()
        # we have k + 1 - p
        for p1 in range(min(k1 + 1, self._p + 1)):
            p2 = self._p - p1
            # print(f"{p1=}, {p2=}, {k1=}, {k2=}")
            # print(f"{vars1=}, {vars2=}")
            # print(f"{new_coeffvars1=} {new_coeffvars2}")
            try:
                ret += self.func(p1, k1, new_genvars1, new_coeffvars1) * self.func(p2, k2, vars1, new_coeffvars2)
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
        return tuple(self.args[2])

    @property
    def coeffvars(self):
        return tuple(self.args[3])

    def _eval_expand_func(self, *args, **kwargs):  # noqa: ARG002
        return elem_sym_poly(self._p, self._k, self.genvars, self.coeffvars)

    @property
    def func(self):
        def e(*args):
            return ElemSym(*args)

        return e

    def xreplace(self, rule):
        res = self
        for v1, v2 in rule.items():
            if v1 in res.genvars:
                new_genvars = [*res.genvars]
                for i, v in enumerate(res.genvars):
                    if v == v1:
                        new_genvars[i] = v2
                res = res.func(res._p, res._k, new_genvars, res.coeffvars)
            if v1 in self.coeffvars:
                new_coeffvars = [*res.coeffvars]
                for i, v in enumerate(res.coeffvars):
                    if v == v1:
                        new_coeffvars[i] = v2
                res = res.func(res._p, res._k, res.genvars, new_coeffvars)
        return res

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
        return printer._print_Function(self)
