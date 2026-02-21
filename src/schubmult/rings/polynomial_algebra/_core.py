from functools import cache

import schubmult.abc as abc
from schubmult.symbolic import (
    EXRAW,
    CoercionFailed,
    S,
    sstr,
    sympify,
)
from schubmult.utils.logging import get_logger

from ..base_ring import BaseRing, BaseRingElement
from ..schubert.base_schubert_ring import BaseSchubertElement
from .polynomial_basis import MonomialBasis

# from .polynomial_basis import EXBasis

logger = get_logger(__name__)


# keys are tuples of nonnegative integers
class PolynomialAlgebraElement(BaseRingElement):

    def __hash__(self):
        return hash(set(self.items()))

    def as_coefficients_dict(self):
        return {self.ring.printing_term(k, self.ring): sympify(v) for k, v in self.items()}

    # def expand(self, deep=True, *args, **kwargs):
    #     return self.change_basis(EXBasis())

    @property
    def free_symbols(self):
        return set()

    def change_basis(self, other_basis):
        new_ring = PolynomialAlgebra(basis=other_basis)
        tfunc = self.ring._basis.transition(other_basis)
        return new_ring.from_dict(tfunc(self))

    def __eq__(self, other):
        if isinstance(other, PolynomialAlgebraElement):
            if other.ring == self.ring:
                diff = self - other
                return all(v == S.Zero for v in diff.values())
        return False

    def expand(self):
        return self.ring._basis.expand(self)

    def __str__(self):
        return sstr(self)


class PolynomialAlgebra(BaseRing):
    def __str__(self):
        return self.__class__.__name__

    def __hash__(self):
        return hash((self.domain, "whatabong"))

    def __eq__(self, other):
        return type(self) is type(other) and self.domain == other.domain and self._basis == other._basis

    def __init__(self, basis, domain=None):
        super().__init__(domain=domain)
        if domain is None:
            self.domain = EXRAW
            self.dom = self.domain

        self._basis = basis
        self.zero_monom = self._basis.zero_monom
        self.dtype = type("PolynomialAlgebraElement", (PolynomialAlgebraElement,), {"ring": self})

    @cache
    def coproduct_on_basis(self, key):
        T = self @ self
        return T.from_dict(self._basis.coproduct(key))

    def _mul_elements(self, elem, other):
        ret = self.zero
        for k0, v0 in elem.items():
            for k, v in other.items():
                ret += self.from_dict(self._basis.product(k0, k, v * v0))
        return ret

    def to_domain(self):
        return self

    def new(self, *x):
        if len(x) == 1:
            if self._basis.is_key(x[0]):
                return self.from_dict({self._basis.as_key(x[0]): S.One})
            return self.from_expr(x[0])
        if self._basis.is_key(x):
            return self.from_dict({self._basis.as_key(x): S.One})
        return self.from_dict(self._basis.from_expr(x))

    def from_expr(self, x):
        return self.from_dict(self._basis.from_expr(x))

    def printing_term(self, k):
        return self._basis.printing_term(k)

    def _coerce_mul(self, other):
        if isinstance(other, PolynomialAlgebraElement) and other.ring == self:
            return other
        return False

    def from_dict(self, element):
        elem = self.dtype()
        elem.update(element)
        return elem

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        if isinstance(element, PolynomialAlgebraElement) or isinstance(element, BaseSchubertElement):
            raise CoercionFailed("Not a domain element")
        return sympify(element)


PA = PolynomialAlgebra(MonomialBasis(abc.x))
