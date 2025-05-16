from functools import cache

from schubmult.rings.poly_lib import complete_sym_poly
from schubmult.rings.variables import NotEnoughGeneratorsError, ZeroGeneratingSet
from schubmult.symbolic import Add, Function, Integer, Mul, Pow, S, sympify, sympify_sympy
from schubmult.symmetric_polynomials.elem_sym import FactorialElemSym
from schubmult.utils.logging import get_logger

from .functions import coeffvars, degree, genvars, is_of_func_type, numvars

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

    def __hash__(self):
        return hash((*self.args, "forp"))

    @property
    def func(self):
        return self.__class__

    @property
    def coeffvars(self):
        return self._coeffvars

    def _sympystr(self, printer):
        return printer._print_Function(self)

    def _eval_expand_func(self, *args, **kwargs):  # noqa: ARG002
        return sympify_sympy(complete_sym_poly(self._p, self._k, self.genvars, self.coeffvars))


class H(CompleteSym_base):
    is_commutative = True
    is_Atom = False
    is_polynomial = True
    is_Function = True
    is_nonzero = True

    @property
    def genvars(self):
        return self._genvars

    @property
    def coeffvars(self):
        return self._coeffvars

    @property
    def degree(self):
        return self._p

    @property
    def numvars(self):
        return self._k

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
        return (S.NegativeOne ** (self._p % 2)) * FactorialElemSym(self._p, self._k + 1 - self._p, self.coeffvars, self.genvars)

    @staticmethod
    def __xnew__(_class, p, k, var1, var2):
        if p < 0 or k < 0:
            return S.Zero
        if p == 0:
            return S.One
        if k == 0:
            return S.Zero
        var1 = var1[:k]
        var2 = var2[: p + k - 1]
        for i, v in enumerate(var1):
            if v in var2:
                j = var2.index(v)
                return FactorialCompleteSym.__new__(_class, p, k - 1, [*var1[:i], *var1[i + 1 :]], [*var2[:j], *var2[j + 1 :]])
        if len(var1) < k:
            raise NotEnoughGeneratorsError(f"{k} passed as number of variables but only {len(var1)} given")
        if len(var2) < p + k - 1:
            raise NotEnoughGeneratorsError(f"{k} passed as number of variables and degree is {p} but only {len(var2)} coefficient variables given. {k + 1 - p} coefficient variables are needed.")
        var1 = tuple(sorted(var1, key=lambda x: sympify_sympy(x).sort_key()))
        var2 = tuple(sorted(var2, key=lambda x: sympify_sympy(x).sort_key()))
        obj = CompleteSym_base.__new__(
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

    @classmethod
    def from_elem_sym(cls, elem, sign=False):
        if sign:
            return (S.NegativeOne ** degree(elem)) * cls(degree(elem), numvars(elem) + 1 - degree(elem), *coeffvars(elem), *genvars(elem))
        return cls(degree(elem), numvars(elem) + 1 - degree(elem), *coeffvars(elem), *genvars(elem))

    @property
    def free_symbols(self):
        return set(self.genvars).union(set(self.coeffvars))

    def split_out_vars(self, vars1, vars2=None):  # noqa: ARG002
        first_vars = [sympify(v) for v in vars1 if v in self.genvars]
        second_vars = [a for a in self.genvars if a not in vars1]
        if len(second_vars) == 0 or len(first_vars) == 0:
            return self
        k1 = len(first_vars)
        k2 = len(second_vars)
        if k1 + k2 != self._k:
            raise Exception

        return Add(
            *[FactorialCompleteSym(i, k1, first_vars, self.coeffvars[: k1 + i - 1]) * FactorialCompleteSym(self._p - i, k2, second_vars, self.coeffvars[k1 + i :]) for i in range(self._p + 1)],
        )

    @property
    def func(self):
        return self.__class__

    def divide_out_diff(self, v1, v2):
        new_obj = self.to_elem_sym().divide_out_diff(v1, v2)
        return new_obj.replace(FactorialElemSym, lambda x: FactorialCompleteSym.from_elem_sym(x))

    @staticmethod
    def from_expr_elem_sym(expr):
        # return expr.replace(FactorialElemSym, lambda *x: (S.NegativeOne**int(x[0]))*FactorialCompleteSym.from_elem_sym(FactorialElemSym(*x)))
        if not expr.args:
            return expr
        if is_of_func_type(expr, FactorialElemSym):
            return FactorialCompleteSym.from_elem_sym(expr, sign=True)
        if isinstance(expr, Mul):
            return Mul(*[FactorialCompleteSym.from_expr_elem_sym(arg) for arg in expr.args])
        if isinstance(expr, Add):
            return Add(*[FactorialCompleteSym.from_expr_elem_sym(arg) for arg in expr.args])
        if isinstance(expr, Pow):
            return Pow(FactorialCompleteSym.from_expr_elem_sym(expr.args[0]), expr.args[1])
        return expr

    @staticmethod
    def to_expr_elem_sym(expr):
        # return expr.replace(FactorialElemSym, lambda *x: (S.NegativeOne**int(x[0]))*FactorialCompleteSym.from_elem_sym(FactorialElemSym(*x)))
        if not expr.args:
            return expr
        if is_of_func_type(expr, FactorialCompleteSym):
            return sympify_sympy(expr).to_elem_sym()
        if isinstance(expr, Mul):
            return Mul(*[FactorialCompleteSym.to_expr_elem_sym(arg) for arg in expr.args])
        if isinstance(expr, Add):
            return Add(*[FactorialCompleteSym.to_expr_elem_sym(arg) for arg in expr.args])
        if isinstance(expr, Pow):
            return Pow(FactorialCompleteSym.to_expr_elem_sym(expr.args[0]), expr.args[1])
        return expr

    def _eval_div_diff(self, v1, v2):
        return FactorialCompleteSym.from_expr_elem_sym(self.to_elem_sym().div_diff(v1, v2))

    def _eval_divide_out_diff(self, v1, v2):
        return FactorialCompleteSym.from_expr_elem_sym(self.to_elem_sym().divide_out_diff(v1, v2))

    def div_diff(self, v1, v2):
        return FactorialCompleteSym.from_expr_elem_sym(self.to_elem_sym().div_diff(v1, v2))


class h(CompleteSym_base):
    is_commutative = True
    is_Atom = False
    is_polynomial = True
    is_Function = True
    is_nonzero = True

    def __new__(cls, p, k, *args):
        p = int(p)
        k = int(k)
        if hasattr(args[0], "__iter__"):
            return CompleteSym.__xnew_cached__(cls, int(p), int(k), tuple(args[0]))
        return CompleteSym.__xnew_cached__(cls, int(p), int(k), tuple(args))

    @staticmethod
    @cache
    def __xnew_cached__(_class, p, k, var1):
        return CompleteSym.__xnew__(_class, p, k, var1)

    @staticmethod
    def __xnew__(_class, p, k, var1):
        if p < 0 or k < 0:
            return S.Zero
        if p == 0:
            return S.One
        if k == 0:
            return S.Zero
        var1 = var1[:k]
        if len(var1) < k:
            raise NotEnoughGeneratorsError(f"{k} passed as number of variables but only {len(var1)} given")
        var1 = tuple(sorted(var1, key=lambda x: sympify_sympy(x).sort_key()))
        obj = CompleteSym_base.__new__(
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

    @property
    def coeffvars(self):
        return ZeroGeneratingSet()


FactorialCompleteSym = H
CompleteSym = h
