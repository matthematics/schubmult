from functools import cache

import schubmult.abc as abc
from schubmult.symbolic import (
    EXRAW,
    Add,
    CoercionFailed,
    CompositeDomain,
    DefaultPrinting,
    DomainElement,
    Ring,
    S,
    sstr,
    sympify,
    sympify_sympy,
    sympy_Add,
    sympy_Mul,
)
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from .base_schubert_ring import BaseSchubertElement
from .polynomial_basis import MonomialBasis

# from .polynomial_basis import EXBasis

logger = get_logger(__name__)


# keys are tuples of nonnegative integers
class PolynomialAlgebraElement(DomainElement, DefaultPrinting, dict):
    precedence = 40

    __sympy__ = True

    def parent(self):
        return self.ring

    def eval(self, *args):
        pass

    def __hash__(self):
        return hash(set(self.items()))

    def _sympystr(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(sympy_Add(*self.as_ordered_terms()), order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def _pretty(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def _latex(self, printer):
        if len(self.keys()) == 0:
            return printer._print(S.Zero)
        if printer.order in ("old", "none"):  # needed to avoid infinite recursion
            return printer._print_Add(self, order="lex")
        return printer._print_Add(sympy_Add(*self.as_ordered_terms()))

    def as_terms(self):
        if len(self.keys()) == 0:
            return [sympify_sympy(S.Zero)]
        return [self[k] if k == () else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k)) for k in self.keys()]

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympify(S.Zero)]
        return [((self[k]) if k == () else sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k))) for k in sorted(self.keys())]

    def __add__(self, other):
        if isinstance(other, PolynomialAlgebraElement):
            if self.ring == other.ring:
                return self.ring.add(self, other)
            return other.__radd__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self._basis.zero_monom: other})
            return self.ring.add(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__add__(new_other)
        except CoercionFailed:
            return other.__radd__(self)

    def __radd__(self, other):
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self._basis.zero_monom: other})
            return self.ring.add(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__add__(self)
        except CoercionFailed:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, PolynomialAlgebraElement):
            if self.ring == other.ring:
                return self.ring.sub(self, other)
            return other.__rsub__(self)
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self._basis.zero_monom: other})
            return self.ring.sub(self, other)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return self.__sub__(new_other)
        except CoercionFailed:
            return other.__rsub__(self)

    def __rsub__(self, other):
        try:
            other = self.ring.domain_new(other)
            other = self.ring.from_dict({self._basis.zero_monom: other})
            return self.ring.sub(other, self)
        except CoercionFailed:
            pass
        try:
            new_other = self.ring(other)
            return new_other.__sub__(self)
        except CoercionFailed:
            return NotImplemented

    def __neg__(self):
        return self.ring.neg(self)

    def __mul__(self, other):
        try:
            return self.ring.mul(self, other)
        except CoercionFailed:
            return other.__rmul__(self)

    def __rmul__(self, other):
        try:
            return self.ring.rmul(self, other)
        except CoercionFailed:
            return NotImplemented

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

    def coproduct(self):
        T = self.ring @ self.ring
        res = T.zero

        for key, val in self.items():
            res += val * self.ring.coproduct_on_basis(key)
        return res

    def as_expr(self):
        return Add(*self.as_terms())

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


class PolynomialAlgebra(Ring, CompositeDomain):
    def __str__(self):
        return self.__class__.__name__

    def __hash__(self):
        return hash((self.domain, "whatabong"))

    def __eq__(self, other):
        return type(self) is type(other) and self.domain == other.domain

    def __init__(self, basis, domain=None):
        if domain:
            self.domain = domain
        else:
            self.domain = EXRAW
        self.dom = self.domain

        self._basis = basis
        self.zero_monom = self._basis.zero_monom
        self.dtype = type("PolynomialAlgebraElement", (PolynomialAlgebraElement,), {"ring": self})

    def __matmul__(self, other):
        from .tensor_ring import TensorRing

        return TensorRing(self, other)

    @cache
    def coproduct_on_basis(self, key):
        T = self @ self
        return T.from_dict(self._basis.coproduct(key))

    def add(self, elem, other):
        return self.from_dict(add_perm_dict(elem, other))

    def sub(self, elem, other):
        return self.from_dict(add_perm_dict(elem, {k: -v for k, v in other.items()}))

    def neg(self, elem):
        return self.from_dict({k: -v for k, v in elem.items()})

    def rmul(self, elem, other):
        if isinstance(other, PolynomialAlgebraElement):
            raise NotImplementedError
        return self.from_dict({k: v * other for k, v in elem.items()})

    def mul(self, elem, other):
        try:
            other = self.domain_new(other)
            return self.from_dict({k: other * v for k, v in elem.items()})
        except Exception:
            pass
        if isinstance(other, PolynomialAlgebraElement):
            ret = self.zero
            for k0, v0 in elem.items():
                for k, v in other.items():
                    ret += self.from_dict(self._basis.product(k0, k, v * v0))
            return ret
        raise CoercionFailed

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

    def _coerce_mul(self, other): ...

    @property
    def one(self):
        return self.from_dict({self._basis.zero_monom: S.One})

    def from_dict(self, element):
        # poly = self.zero
        # for monom, coeff in element.items():
        #     if coeff != self.domain.zero:
        #         poly[monom] = coeff
        # return poly
        return self.dtype(element)

    @property
    def zero(self):
        return self.dtype()

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        if isinstance(element, PolynomialAlgebraElement) or isinstance(element, BaseSchubertElement):
            raise CoercionFailed("Not a domain element")
        return sympify(element)


PA = PolynomialAlgebra(MonomialBasis(abc.x))
