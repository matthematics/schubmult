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
    """Element of a polynomial algebra, stored as a dict mapping basis keys to coefficients.

    Keys are exponent tuples (in the monomial basis) or basis-specific keys
    depending on the parent ring's basis. Supports arithmetic, basis changes,
    and duality pairing with free algebra elements.
    """
    def __hash__(self):
        return hash(set(self.items()))

    def as_coefficients_dict(self):
        """Return a dict mapping printing terms to sympified coefficients."""
        return {self.ring.printing_term(k, self.ring): sympify(v) for k, v in self.items()}

    # def expand(self, deep=True, *args, **kwargs):
    #     return self.change_basis(EXBasis())

    @property
    def free_symbols(self):
        return set()

    def change_basis(self, other_basis: type):
        """Convert this element to another polynomial basis.

        Args:
            other_basis: A basis class, basis instance, or callable returning a basis.

        Returns:
            A new PolynomialAlgebraElement in the target basis's ring.
        """
        if isinstance(other_basis, type):
            basis_obj = other_basis(self.ring.genset)
        elif hasattr(other_basis, "transition") and hasattr(other_basis, "is_key"):
            basis_obj = other_basis
        else:
            basis_obj = other_basis(self.ring.genset)
        new_ring = PolynomialAlgebra(basis=basis_obj)
        tfunc = self.ring._basis.transition(new_ring._basis)
        return new_ring.from_dict(tfunc(self))

    def __eq__(self, other):
        if isinstance(other, PolynomialAlgebraElement):
            if other.ring == self.ring:
                diff = self - other
                return all(v == S.Zero for v in diff.values())
        return False

    def expand(self):
        """Expand this element into an explicit polynomial expression."""
        return self.ring._basis.expand(self)

    def apply_dual_element(self, dual_elem):
        """Pair this polynomial element with a dual free algebra element.

        Converts *self* to the monomial basis and *dual_elem* to the word
        basis, then sums products of matching coefficients.

        Args:
            dual_elem: A FreeAlgebraElement to pair with.

        Returns:
            The scalar pairing value.
        """
        from ..free_algebra.word_basis import WordBasis
        from .monomial_basis import MonomialBasis

        word_elem = dual_elem.change_basis(WordBasis)
        monom_elem = self.change_basis(MonomialBasis(self.ring.genset))

        result = S.Zero

        for k, v in monom_elem.items():
            for k2, v2 in word_elem.items():
                if k == k2:
                    result += v * v2
        return result

    def __str__(self):
        return sstr(self)


class PolynomialAlgebra(BaseRing):
    """Polynomial algebra ring with a configurable basis.

    The algebra operates on :class:`PolynomialAlgebraElement` instances whose
    keys are determined by the chosen basis. Supports multiplication, basis
    changes, coproducts, and conversion from symbolic expressions.

    Args:
        basis: A basis instance (e.g. ``MonomialBasis(x)``).
        domain: Coefficient domain (default ``EXRAW``).
    """
    def __str__(self):
        return self.__class__.__name__

    def __hash__(self):
        return hash((self.domain, "whatabong"))

    def __eq__(self, other):
        return type(self) is type(other) and self.domain == other.domain and self._basis == other._basis

    def __init__(self, basis, domain=None):
        """Initialize a PolynomialAlgebra with the given basis and coefficient domain."""
        super().__init__(domain=domain)
        if domain is None:
            self.domain = EXRAW
            self.dom = self.domain

        self._basis = basis
        self.zero_monom = self._basis.zero_monom
        self.dtype = type("PolynomialAlgebraElement", (PolynomialAlgebraElement,), {"ring": self})

    @property
    def genset(self):
        return self._basis.genset

    @cache
    def coproduct_on_basis(self, key):
        """Compute the coproduct of a single basis key in the tensor ring."""
        T = self @ self
        return T.from_dict(self._basis.coproduct(key))

    def mul(self, elem, other):
        """Multiply two elements via the basis product rule."""
        try:
            ret = self.zero
            for k0, v0 in elem.items():
                for k, v in other.items():
                    ret += self.from_dict(self._basis.product(k0, k, v * v0))
            return ret
        except Exception:
            return super().mul(elem, other)

    def to_domain(self):
        return self

    def new(self, *x):
        """Create a new element from the given key or expression."""
        if len(x) == 1:
            if isinstance(x[0], int):
                return self.from_dict({x: S.One})
            if self._basis.is_key(x[0]):
                return self.from_dict({self._basis.as_key(x[0]): S.One})
            return self.from_expr(x[0])
        if self._basis.is_key(x):
            return self.from_dict({self._basis.as_key(x): S.One})
        return self.from_expr(x)

    def from_expr(self, x, length=None):
        """Create an element from a symbolic expression.

        Parses *x* into monomials, then transitions to this ring's basis.

        Args:
            x: A symbolic polynomial expression.
            length: Optional fixed number of variables.

        Returns:
            A PolynomialAlgebraElement in this ring.
        """
        from .monomial_basis import MonomialBasis

        monomial_basis = MonomialBasis(genset=self.genset)
        monomial_dict = monomial_basis.from_expr(x, length=length)

        if isinstance(self._basis, MonomialBasis):
            return self.from_dict(monomial_dict)

        return self.from_dict(monomial_basis.transition(self._basis)(monomial_dict))

    def printing_term(self, k):
        """Return the display symbol for basis key *k*."""
        return self._basis.printing_term(k)

    def _coerce_mul(self, other):
        if isinstance(other, PolynomialAlgebraElement) and other.ring == self:
            return other
        return False

    def from_dict(self, element):
        """Construct an element from a dict of ``{key: coefficient}`` pairs."""
        elem = self.dtype()
        elem.update(element)
        return elem

    def domain_new(self, element, orig_domain=None):  # noqa: ARG002
        """Coerce a raw value into the coefficient domain."""
        if isinstance(element, PolynomialAlgebraElement) or isinstance(element, BaseSchubertElement):
            raise CoercionFailed("Not a domain element")
        return sympify(element)


PA = PolynomialAlgebra(MonomialBasis(abc.x))
