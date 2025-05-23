from functools import cache

from schubmult.rings.poly_lib import elem_sym_poly_q
from schubmult.rings.variables import GeneratingSet, NotEnoughGeneratorsError
from schubmult.symbolic import Integer, S, sympify_sympy

from .elem_sym import ElemSym_base


class E_q(ElemSym_base):
    def __new__(cls, p, k, *args, x_var= None, q_var=None):
        p = int(p)
        k = int(k)
        if not q_var:
            q_var = GeneratingSet("q")
        if not x_var:
            x_var = GeneratingSet("x")
        if hasattr(args[0], "__iter__"):
            return QFactorialElemSym.__xnew_cached__(cls, int(p), int(k), tuple(args[0]), tuple(args[1]), x_var=x_var, q_var=q_var)
        return QFactorialElemSym.__xnew_cached__(cls, int(p), int(k), tuple(args[: int(k)]), tuple(args[k : 2 * k + 1 - p]), x_var = x_var, q_var=q_var)

    def __hash__(self):
        return hash((self._p, self._k, self._genvars, self._coeffvars, self._x_var, self._q_var))

    @property
    def free_symbols(self):
        return {*self._genvars, *self._coeffvars}

    @staticmethod
    @cache
    def __xnew_cached__(_class, p, k, var1, var2, x_var, q_var):
        return QFactorialElemSym.__xnew__(_class, p, k, var1, var2, x_var, q_var)

    @staticmethod
    def __xnew__(_class, p, k, var1, var2, x_var, q_var):
        if p > k or k < 0:
            return S.Zero
        if p == 0:
            return S.One
        var1 = var1[:k]
        var2 = var2[: k + 1 - p]
        for i, v in enumerate(var1):
            if v in var2:
                j = var2.index(v)
                return QFactorialElemSym.__new__(_class, p, k - 1, [*var1[:i], *var1[i + 1 :]], [*var2[:j], *var2[j + 1 :]], q_var=q_var)
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
        obj._q_var = q_var
        obj._x_var = x_var
        return obj

    @cache
    def _eval_expand_func(self, *args, **_):  # noqa: ARG002
        return sympify_sympy(elem_sym_poly_q(self._p, self._k, self.genvars, self.coeffvars, q_var=(0, *[self._q_var[self._x_var.index(b)] for b in self.genvars])))

    # def pull_out_vars(self, var1, var2, min_degree=1):
    #     if self._p < min_degree:
    #         return self
    #     if var1 in self.genvars and var2 in self.coeffvars:
    #         return self.xreplace({var1: var2}) + QFactorialElemSym(1, 1, [var1], [var2]) * self.divide_out_diff(var1, var2)
    #     return self


QFactorialElemSym = E_q
