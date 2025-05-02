from functools import cache

from sympy import Add, Dict, Integer, S, Tuple, Wild, sympify
from sympy.core.expr import Expr

from schubmult.utils.logging import get_logger
from schubmult.utils.ring_utils import NotEnoughGeneratorsError

from .poly_lib import elem_sym_poly

logger = get_logger(__name__)


class ElemSym(Expr):
    is_commutative = True
    is_Atom = False
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
        obj = Expr.__new__(
            _class,
            Integer(p),
            Integer(k),
            Tuple(*sorted(var1, key=lambda x: int(x.name.split("_")[1]))),
            Tuple(*sorted(var2, key=lambda x: int(x.name.split("_")[1]))),
        )
        obj._p = p
        obj._k = k
        obj._genvars = var1
        obj._coeffvars = var2

        return obj

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
        return tuple(self.args[2])

    @property
    def coeffvars(self):
        return tuple(self.args[3])

    def _eval_expand_func(self, *args, **kwargs):  # noqa: ARG002
        return sympify(elem_sym_poly(self._p, self._k, self.genvars, self.coeffvars))

    @property
    def func(self, *args):
        def e(*args):
            return self.__class__(*args)

        return e

    #     return e
    def _eval_subs(self, *rule):
        # print(f"_eval_subs")
        # print(f"{rule=}")
        # print(f"{self=}")
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
        # print(f"xreplace")
        # print(f"{rule=}")
        # print(f"{self=}")
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

    # def __mul__(self, other):
    #     if self._p == 1:
    #         if isinstance(other, ElemSym):
    #             for i, v in enumerate(self.coeffvars):
    #                 if v in other.coeffvars:
    #                     j = 0
    #                     while j < self._k and self.genvars[j] in other.genvars:
    #                         j += 1
    #                     if j == self._k:
    #                         j -= 1
    #                     v1 = self.genvars[j]
    #                     newgenvars1 = [*self.genvars[:j], *self.genvars[j+1:]]
    #                     newcoeffvars1 = [*self.coeffvars[:i],*self.coeffvars[i+1:]]
    #                     newcoeffvars2 = [*other.coeffvars]
    #                     newcoeffvars2.remove(v)
    #                     if self._k == 1:
    #                         return ElemSym(other._p + 1, other._k, other.genvars, newcoeffvars2) - ElemSym(other._p + 1, other._k + 1, [*other.genvars, self.genvars[0]], other.coeffvars)
    #                     return ElemSym(1, self._k - 1, newgenvars1, newcoeffvars1)* (ElemSym(other._p + 1, other._k, other.genvars, newcoeffvars2) - ElemSym(other._p + 1, other._k + 1, [*other.genvars, self.genvars[0]], other.coeffvars))
    #     return Mul(self, other)

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
        return printer._print_Function(self)

    def pull_out_vars(self, var1, var2, min_degree=1):
        if self._p < min_degree:
            return self
        if var1 in self.genvars and var2 in self.coeffvars:
            return self.xreplace({var1: var2}) + ElemSym(1, 1, [var1], [var2]) * self.divide_out_diff(var1, var2)
        return self


def split_out_vars(expr, vars1, vars2):
    expr = sympify(expr)
    if not expr.args:
        return expr
    if hasattr(expr, "split_out_vars"):
        try:
            return expr.split_out_vars(vars1, vars2)
        except NotEnoughGeneratorsError:
            return expr
    return expr.func(*[split_out_vars(arg, vars1, vars2) for arg in expr.args])


def pull_out_vars(expr, var1, var2, min_degree=1):
    expr = sympify(expr)
    if hasattr(expr, "pull_out_vars"):
        return expr.pull_out_vars(var1, var2, min_degree)
    if not expr.args:
        return expr
    return expr.func(*[pull_out_vars(arg, var1, var2, min_degree) for arg in expr.args])
    # return expr.xreplace({sympify(var1): sympify(var2)}) + ElemSym(1, 1, [var1],[var2]) * divide_out_diff(expr, sympify(var1), sympify(var2))


#
# ElemSym(p, k, stuff1 + v1, stuff1 + v2) = ElemSym(p, k, stuff1, stuff2) + Elem(1, 1, v1, v2)*ElemSym(p - 1, k - 1, stuff1, stuff2 + v2)

# def elem_sym_unify_recurse(expr, expr1)


def to_complete_sym(self):
    return S.NegativeOne ** (self._p % 2) * CompleteSym.from_elem_sym(self)


def elem_sym_unify(expr, arg=None):
    expr = sympify(expr)

    if not expr.args:
        return expr
    if isinstance(expr, ElemSym):
        return expr

    if not arg:
        for arg in expr.args:
            expr = elem_sym_unify(expr, arg)
        return expr.func(*[elem_sym_unify(arg) for arg in expr.args])
    expr2 = expr
    if isinstance(arg, ElemSym):
        v1 = Wild("v_1")
        v2 = Wild("v_2")
        rep_pattern = ElemSym(arg._p, arg._k + 1, [*arg.genvars, v1], [*arg.coeffvars, v2])
        pattern = arg + ElemSym(1, 1, [v1], [v2]) * ElemSym(arg._p - 1, arg._k, arg.genvars, [*arg.coeffvars, v2])
        expr2 = expr.replace(pattern, rep_pattern)
    if arg.args:
        for arg2 in arg.args:
            expr2 = elem_sym_unify(expr2, arg2)
    if expr != expr2:
        return elem_sym_unify(expr2)
    return expr


class CompleteSym(Expr):
    is_commutative = True
    is_Atom = False
    is_polynomial = True
    is_Function = True
    is_nonzero = True

    def __new__(cls, p, k, var1, var2):
        return CompleteSym.__xnew_cached__(cls, p, k, tuple(var1), tuple(var2))

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
        obj = Expr.__new__(_class, cont_obj.args[0], cont_obj.args[1] + 1 - cont_obj.args[0], cont_obj.args[3], cont_obj.args[2])
        obj._under_elem = cont_obj
        obj._p = obj.args[0]
        obj._k = obj.args[1]
        return obj

    @classmethod
    def from_elem_sym(cls, elem):
        return cls(elem._p, elem._k + 1 - elem._p, elem.coeffvars, elem.genvars)

    @property
    def free_symbols(self):
        return set(self.genvars).union(set(self.coeffvars))

    def split_out_vars(self, vars1, vars2=None):
        if not all(v in self.genvars for v in vars1):
            return self
        first_vars = [*vars1]
        second_vars = [a for a in self.genvars if a not in vars1]
        k1 = len(first_vars)
        k2 = len(second_vars)
        if k1 + k2 != self._k:
            raise Exception

        return Add(*[CompleteSym(i, k1, first_vars, self.coeffvars[: k1 + i - 1]) * CompleteSym(self._p - i, k2, second_vars, self.coeffvars[k1 + i :]) for i in range(0, self._p + 1)])

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

    def _eval_expand_func(self, *args, **kwargs):
        return sympify(elem_sym_poly(self._under_elem._p, self._under_elem._k, [-x for x in self._under_elem.genvars], [-y for y in self._under_elem.coeffvars]))

    @property
    def func(self, *args):
        def h(*args):
            return self.__class__(*args)

        return h

    #     return e
    def _eval_subs(self, *rule):
        # print(f"_eval_subs")
        # print(f"{rule=}")
        # print(f"{self=}")
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
        # print(f"xreplace")
        # print(f"{rule=}")
        # print(f"{self=}")
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
    def from_expr_elem_sym(self, expr):
        return expr.replace(ElemSym, lambda x: CompleteSym.from_elem_sym(ElemSym(*x)))

    def _eval_div_diff(self, v1, v2):
        return CompleteSym.from_expr_elem_sym(self._under_elem.div_diff(v1, v2))

    def _eval_divide_out_diff(self, v1, v2):
        return CompleteSym.from_expr_elem_sym(self._under_elem.divide_out_diff(v1, v2))

    def div_diff(self, v1, v2):
        return CompleteSym.from_expr_elem_sym(self._under_elem.div_diff(v1, v2))

    def _sympystr(self, printer):
        return printer._print_Function(self)
