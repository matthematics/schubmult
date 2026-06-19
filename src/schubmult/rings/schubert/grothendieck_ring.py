from functools import cache

import schubmult.rings.printing as spolymod
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, Symbol
from schubmult.symbolic.poly.schub_poly import groth_mul_full, grothendieck_poly
from schubmult.symbolic.poly.variables import GeneratingSet, ZeroGeneratingSet

from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing

__all__ = [
    "GrothendieckElement",
    "GrothendieckRing",
    "Gx",
]


class GrothendieckElement(BaseSchubertElement):
    """Element of a GrothendieckRing, stored as {Permutation: coeff}."""

    def as_polynomial(self):
        from schubmult.symbolic import Add
        return Add(*[v * self.ring.cached_schubpoly(k) for k, v in self.items()])


class GrothendieckRing(BaseSchubertRing):
    """
    Ring of Grothendieck polynomials.

    A deformation of the Schubert polynomial ring with parameter beta.
    Basis elements G_w satisfy G_u * G_v = sum_w c^w_{u,v}(beta) G_w
    where c^w_{u,v}(beta) are polynomials in beta with integer coefficients.
    When beta=0, recovers ordinary Schubert polynomials.

    Parameters
    ----------
    genset : GeneratingSet
        The generating set (variable alphabet).
    beta : sympy/symengine symbol, optional
        The deformation parameter. Defaults to Symbol("β").
    """

    def __init__(self, genset, beta=None):
        super().__init__(genset, coeff_genset=None)
        if beta is None:
            beta = Symbol("\u03B2")
        self._beta = beta
        self._zz = ZeroGeneratingSet()
        # Rebind dtype to GrothendieckElement subclass
        self.dtype = type("GrothendieckElement", (GrothendieckElement,), {"ring": self})

    def __str__(self):
        return f"Grothendieck polynomial ring in {self.genset.label}"

    def __repr__(self):
        return f"GrothendieckRing({self.genset!r}, beta={self._beta!r})"

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "GR", self._beta))

    def __eq__(self, other):
        return type(self) is type(other) and self.genset == other.genset and self._beta == other._beta

    @property
    def beta(self):
        return self._beta

    @cache
    def cached_product(self, u, v, basis2):
        if self == basis2:
            return groth_mul_full({u: S.One}, v, self.genset, self._zz, self._beta)
        # Cross-ring fallback: delegate to single Schubert multiplication
        return super().cached_product(u, v, basis2)

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    @cache
    def cached_schubpoly(self, k):
        return grothendieck_poly(k, self.genset, self._zz, self._beta)

    def printing_term(self, k, prefix=""):
        return spolymod.GrothendieckPoly(k, self.genset.label, prefix=prefix)

    def __call__(self, x):
        return self.new(x)

    def new(self, x):
        genset = self.genset
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self.from_dict({Permutation(x): self.domain.one})
        elif isinstance(x, Permutation):
            elem = self.from_dict({x: S.One})
        elif isinstance(x, GrothendieckElement):
            if x.ring.genset == genset:
                return x
            raise ValueError("Different generating set")
        else:
            elem = self.from_expr(x)
        return elem

    def from_dict(self, dct):
        return self.dtype(dct)


Gx = GrothendieckRing(GeneratingSet("x"))
