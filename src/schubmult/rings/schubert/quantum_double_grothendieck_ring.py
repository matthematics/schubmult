"""Quantum double Grothendieck polynomials: the ``QDGx`` interface.

Multiplication uses the conjectural equivariant quantum K-Pieri kernel;
explicit polynomial expansion uses Lenart--Maeno quantization. This is the
stable polynomial ring: terms outside a chosen finite symmetric group are
not discarded.
"""

from functools import cache

from schubmult.abc import beta as default_beta
from schubmult.combinatorics.permutation import Permutation
from schubmult.mult.groth_double import normalize_coeff
from schubmult.mult.groth_quantum_double import grothmult_q_double, qgroth_poly
from schubmult.rings.base_ring import BaseRingElement
from schubmult.rings.printing import QuantumDoubleGrothendieckPoly
from schubmult.rings.schubert.base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from schubmult.symbolic import CoercionFailed, S, expand, sympify, sympify_sympy, sympy_Add
from schubmult.symbolic.poly.variables import GeneratingSet, ZeroGeneratingSet

__all__ = ["QDGx", "QuantumDoubleGrothendieckElement", "QuantumDoubleGrothendieckRing"]

_zero_coeff_genset = ZeroGeneratingSet()


class QuantumDoubleGrothendieckElement(BaseSchubertElement):
    """A linear combination of quantum double Grothendieck polynomials."""

    def __eq__(self, other):
        return isinstance(other, QuantumDoubleGrothendieckElement) and self.ring == other.ring and dict.__eq__(self, other)

    def __ne__(self, other):
        return not self == other

    def as_expr(self):
        """Return display atoms through SymPy, which supports the ring's printing terms."""
        return sympy_Add(*self.as_terms())

    def simplify(self, factor=True):
        """Cancel coefficients, optionally factoring their numerators."""
        coefficients = {}
        for perm, coefficient in self.items():
            value = expand(coefficient) if self.ring.beta == 0 else normalize_coeff(coefficient, self.ring.coeff_genset, beta=self.ring.beta)
            if value == S.Zero:
                continue
            if factor:
                numerator, denominator = sympify_sympy(value).as_numer_denom()
                value = sympify(numerator.factor() / denominator)
            coefficients[perm] = value
        return self.ring.from_dict(coefficients)

    def mult_poly(self, poly):
        """Multiply by a coefficient; conversion of arbitrary x-polynomials is not implemented."""
        return self.ring.mul_expr(self, poly)


class QuantumDoubleGrothendieckRing(BaseSchubertRing):
    """Stable quantum double Grothendieck polynomials with parameter ``beta``.

    Products, including mixed coefficient alphabets, use ``grothmult_q_double``
    directly. ``as_polynomial()`` and ``expand()`` use ``qgroth_poly``.
    Construction from arbitrary expressions involving x is not implemented;
    construct basis elements from permutations instead.
    """

    def __init__(self, genset, coeff_genset, beta=None):
        if isinstance(coeff_genset, ZeroGeneratingSet):
            coeff_genset = _zero_coeff_genset
        super().__init__(genset, coeff_genset)
        self._beta = default_beta if beta is None else sympify(beta)
        self.dtype = type("QuantumDoubleGrothendieckElement", (QuantumDoubleGrothendieckElement,), {"ring": self})

    @property
    def beta(self):
        """The K-theoretic deformation parameter."""
        return self._beta

    def __str__(self):
        return f"Quantum double Grothendieck polynomial ring in {self.genset.label} and {self.coeff_genset.label}"

    def __repr__(self):
        return f"QuantumDoubleGrothendieckRing({self.genset!r}, {self.coeff_genset!r}, beta={self.beta!r})"

    def __eq__(self, other):
        return type(self) is type(other) and self.genset == other.genset and self.coeff_genset == other.coeff_genset and self.beta == other.beta

    def __hash__(self):
        return hash((type(self), self.genset, self.coeff_genset, self.beta))

    def __call__(self, value):
        return self.new(value)

    def new(self, value):
        """Construct from a permutation, an element of this ring, or a coefficient."""
        if isinstance(value, Permutation | list | tuple):
            return self.from_dict({Permutation(value): S.One})
        if isinstance(value, BaseRingElement):
            if value.ring == self:
                return value
            raise ValueError(f"Cannot convert an element of {value.ring} to {self}")
        return self.from_expr(value)

    def from_dict(self, coefficients):
        """Construct from permutation coefficients, dropping exact zeros."""
        return self.dtype({Permutation(w): sympify(c) for w, c in coefficients.items() if c != S.Zero})

    def add(self, elem, other):
        if other.ring != self:
            raise ValueError("Addition requires identical quantum double Grothendieck rings")
        return super().add(elem, other)

    def sub(self, elem, other):
        if other.ring != self:
            raise ValueError("Subtraction requires identical quantum double Grothendieck rings")
        return super().sub(elem, other)

    def mul(self, elem, other):
        """Multiply through the quantum double Grothendieck kernel."""
        if isinstance(other, BaseRingElement):
            if not isinstance(other, QuantumDoubleGrothendieckElement) or self.genset != other.ring.genset or self.beta != other.ring.beta:
                raise ValueError("Multiplication requires quantum double Grothendieck rings with identical x alphabets and beta")
        return super().mul(elem, other)

    @cache
    def cached_product(self, u, v, basis2):
        """Kernel structure constants in the left factor's coefficient alphabet."""
        if not isinstance(basis2, QuantumDoubleGrothendieckRing) or self.genset != basis2.genset or self.beta != basis2.beta:
            raise ValueError("Incompatible quantum double Grothendieck rings")
        return grothmult_q_double({u: S.One}, v, self.coeff_genset, basis2.coeff_genset, beta=self.beta)

    @cache
    def cached_positive_product(self, u, v, basis2):
        """Use the same kernel, without a separate positivity transformation."""
        return self.cached_product(u, v, basis2)

    @cache
    def cached_schubpoly(self, perm):
        """The explicit Lenart--Maeno quantum double Grothendieck polynomial."""
        # The existing quantization divides by beta before cancelling.
        beta = default_beta if self.beta == 0 else self.beta
        polynomial = qgroth_poly(perm, self.genset, self.coeff_genset, beta=beta)
        return polynomial.subs(default_beta, 0) if self.beta == 0 else polynomial

    def printing_term(self, perm, prefix=""):
        return QuantumDoubleGrothendieckPoly(perm, self.genset.label, self.coeff_genset.label, prefix=f"{prefix}Q")

    def mul_expr(self, elem, expr):
        """Multiply by a coefficient, refusing unsupported x-polynomial conversion."""
        try:
            coefficient = self.domain_new(expr)
        except CoercionFailed as error:
            raise NotImplementedError("Conversion of x-polynomials to the quantum Grothendieck basis is not implemented") from error
        return self.from_dict({w: coefficient * c for w, c in elem.items()})

    def in_quantum_basis(self, elem):
        return elem

    def in_classical_basis(self, elem):
        """Re-express the explicit polynomial in the classical double Grothendieck basis."""
        from schubmult.rings.schubert.double_grothendieck_ring import DoubleGrothendieckRing

        return DoubleGrothendieckRing(self.genset, self.coeff_genset, beta=self.beta).from_expr(elem.as_polynomial())


def QDGx(value, genset=GeneratingSet("y"), beta=None):
    """Construct a quantum double Grothendieck element; ``genset`` selects y (or z).

    ``beta=-1`` uses the usual K-theory convention; ``genset="0"`` gives
    non-equivariant quantum Grothendieck polynomials.
    """
    if isinstance(genset, str):
        genset = ZeroGeneratingSet() if genset == "0" else GeneratingSet(genset)
    return QuantumDoubleGrothendieckRing(GeneratingSet("x"), genset, beta=beta)(value)
