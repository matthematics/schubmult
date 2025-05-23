from functools import cache

from schubmult.rings.poly_lib import elem_sym_poly
from schubmult.rings.variables import NotEnoughGeneratorsError, ZeroGeneratingSet
from schubmult.symbolic import Add, Function, Integer, S, sympify, sympify_sympy
from schubmult.utils.logging import get_logger

logger = get_logger(__name__)

# class _E(SympyExpr):
#     pass

# class Buns:
#     pass

# class Sponge(SymengineExprClass):
#     @property
#     def __class__(self):
#         return Buns


class ElemSym_base(Function):
    is_commutative = True
    is_Atom = False
    is_polynomial = True
    is_Function = True
    is_nonzero = True

    def _sympystr(self, printer):
        return printer._print_Function(self)

    def _eval_subs(self, old, new):
        rule = {old: new}
        new_args = [*self.args]
        for i, arg in enumerate(self.args[2:]):
            new_args[i + 2] = sympify(arg).subs(rule)
        return self.func(*new_args)

    def subs(self, *rule):
        rule = dict(rule)
        new_args = [*self.args]
        for i, arg in enumerate(self.args[2:]):
            new_args[i + 2] = sympify(arg).subs(rule)
        return self.func(*new_args)

    def xreplace(self, rule):
        new_args = [*self.args]
        for i, arg in enumerate(self.args[2:]):
            new_args[i + 2] = sympify(arg).xreplace(rule)
        return self.func(*new_args)

    def has(self, *args):
        return any(arg in self.genvars for arg in args) or any(arg in self.coeffvars for arg in args)

    def has_symbol(self, sym):
        return sym in self.genvars or sym in self.coeffvars

    @property
    def is_Symbol(self):
        return False

    @property
    def is_symbol(self):
        return False

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


    @property
    def func(self):
        return self.__class__


class E(ElemSym_base):
    def __new__(cls, p, k, *args):
        p = int(p)
        k = int(k)
        if hasattr(args[0], "__iter__"):
            return FactorialElemSym.__xnew_cached__(cls, int(p), int(k), tuple(args[0]), tuple(args[1]))
        return FactorialElemSym.__xnew_cached__(cls, int(p), int(k), tuple(args[: int(k)]), tuple(args[k : 2 * k + 1 - p]))

    @property
    def free_symbols(self):
        return {*self._genvars, *self._coeffvars}

    @staticmethod
    @cache
    def __xnew_cached__(_class, p, k, var1, var2):
        return FactorialElemSym.__xnew__(_class, p, k, var1, var2)

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
                return FactorialElemSym.__new__(_class, p, k - 1, [*var1[:i], *var1[i + 1 :]], [*var2[:j], *var2[j + 1 :]])
        if len(var1) < k:
            raise NotEnoughGeneratorsError(f"{k} passed as number of variables but only {len(var1)} given")
        if len(var2) < k + 1 - p:
            raise NotEnoughGeneratorsError(f"{k} passed as number of variables and degree is {p} but only {len(var2)} coefficient variables given. {k + 1 - p} coefficient variables are needed.")
        var1 = tuple(sorted(var1, key=lambda x: sympify_sympy(x).sort_key()))
        var2 = tuple(sorted(var2, key=lambda x: sympify_sympy(x).sort_key()))
        obj = ElemSym_base.__new__(
            _class,
            Integer(p),
            Integer(k),
            *var1,
            *var2,
        )
        obj._p = p
        obj._k = k
        obj._genvars = var1
        obj._coeffvars = var2
        return obj

    def split_out_vars(self, vars1, vars2=None):
        if vars1 is None:
            first_vars = [sympify(v) for v in vars2 if v in self.coeffvars]
            second_vars = [a for a in self.coeffvars if a not in vars2]
            k1 = len(first_vars)
            k2 = len(second_vars)
            # if k1 + k2 != self._k + 1 - self.:
            #     raise Exception

            return Add(
                *[
                    FactorialElemSym(i, k1 + i - 1, self.genvars[: k1 + i - 1], first_vars) * FactorialElemSym(self._p - i, len(self.genvars[k1 + i :]), self.genvars[k1 + i :], second_vars)
                    for i in range(self._p + 1)
                ],
            )
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

        ret = S.Zero
        k1 = len(vars1)
        k2 = self._k - k1
        new_genvars1 = [*self.genvars]
        for v in vars1:
            new_genvars1.remove(v)
        new_coeffvars2 = [*newvars2]
        new_coeffvars1 = [*newvars2]
        new_coeffvars2.reverse()
        for p1 in range(min(k1 + 1, self._p + 1)):
            p2 = self._p - p1
            try:
                ret += self.func(p1, k1, vars1, new_coeffvars1) * self.func(p2, k2, new_genvars1, new_coeffvars2)
            except NotEnoughGeneratorsError:
                pass
        return ret

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

    def pull_out_vars(self, var1, var2, min_degree=1):
        if self._p < min_degree:
            return self
        if var1 in self.genvars and var2 in self.coeffvars:
            return self.xreplace({var1: var2}) + FactorialElemSym(1, 1, [var1], [var2]) * self.divide_out_diff(var1, var2)
        return self

    @cache
    def _eval_expand_func(self, *args, **_):  # noqa: ARG002
        return sympify_sympy(elem_sym_poly(self._p, self._k, self.genvars, self.coeffvars))


class e(ElemSym_base):
    def __new__(cls, p, k, *args):
        p = int(p)
        k = int(k)
        if hasattr(args[0], "__iter__"):
            return ElemSym.__xnew_cached__(cls, p, k, tuple(args[0]))
        return ElemSym.__xnew_cached__(cls, p, k, tuple(args))

    @staticmethod
    @cache
    def __xnew_cached__(_class, p, k, var1):
        return ElemSym.__xnew__(_class, p, k, var1)

    @staticmethod
    def __xnew__(_class, p, k, var1):
        if p > k or k < 0:
            return S.Zero
        if p == 0:
            return S.One
        var1 = var1[:k]
        if len(var1) < k:
            raise NotEnoughGeneratorsError(f"{k} passed as number of variables but only {len(var1)} given")
        var1 = tuple(sorted(var1, key=lambda x: sympify_sympy(x).sort_key()))
        obj = ElemSym_base.__new__(
            _class,
            Integer(p),
            Integer(k),
            *var1,
        )
        obj._p = p
        obj._k = k
        obj._genvars = var1

        return obj

    @property
    def free_symbols(self):
        return set(self.genvars)

    def split_out_vars(self, vars1, vars2=None):
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
        for p1 in range(min(k1 + 1, self._p + 1)):
            p2 = self._p - p1
            try:
                ret += self.func(p1, k1, vars1) * self.func(p2, k2, new_genvars1)
            except NotEnoughGeneratorsError:
                pass
        return ret

    @property
    def coeffvars(self):
        return ZeroGeneratingSet()

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

    @cache
    def _eval_expand_func(self, *args, **_):  # noqa: ARG002
        return sympify_sympy(elem_sym_poly(self._p, self._k, self.genvars, self.coeffvars))


FactorialElemSym = E

ElemSym = e
