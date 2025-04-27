import sympy
from symengine import expand, sympify


# Atomic Schubert polynomial
class AbstractSchubPoly(sympy.Expr):
    is_Atom = True

    def __new__(cls, k, basis):
        obj = sympy.Expr.__new__(cls, k, basis)
        obj._key = k
        obj._genset = basis.genset
        obj._basis = basis
        obj._perm = k
        return obj

    def _eval_expand_basic(self, *args, **kwargs):  # noqa: ARG002
        return self.as_polynomial()

    def expand(self, deep=True, *args, **kwargs):  # noqa: ARG002
        if not deep:
            return self.ring.from_dict({k: expand(v) for k, v in self.items()})
        return sympy.sympify(expand(sympify(self.as_polynomial())))

    def as_polynomial(self):
        return self.ring.cached_schubpoly(self._key)

    @property
    def ring(self):
        return self._basis

    @property
    def genset(self):
        return self._genset

    @property
    def coeff_genset(self):
        return self.ring.coeff_genset

    @property
    def args(self):
        return (sympy.Tuple(*self._key), self._basis)
