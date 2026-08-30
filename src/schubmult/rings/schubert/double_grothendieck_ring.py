from functools import cache, cached_property

import schubmult.rings.printing as spolymod
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, Symbol, sympify
from schubmult.symbolic.common_polys import grothendieck_poly_with_ring
from schubmult.symbolic.poly.variables import GeneratingSet

from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from .double_schubert_ring import DoubleSchubertRing

__all__ = [
    "DGx",
    "DoubleGrothendieckElement",
    "DoubleGrothendieckRing",
]


class DoubleGrothendieckElement(BaseSchubertElement):
    """Element of a DoubleGrothendieckRing, stored as {Permutation: coeff}."""

    def as_polynomial(self):
        from schubmult.symbolic import Add

        return Add(*[v * self.ring.cached_schubpoly(k) for k, v in self.items()])

    def perm_subs(self, perm):
        return self.ring.perm_subs(self, perm)


class DoubleGrothendieckRing(BaseSchubertRing):
    """
    Ring of double (K-theoretic) Grothendieck polynomials G_w(x, y).

    Elements are stored in the G-basis as ``{Permutation: coeff}``. There is
    no direct structure-constant formula for the product implemented here:
    instead both factors are expanded into the underlying ``DoubleSchubertRing``
    (via ``grothendieck_poly_with_ring``), multiplied there, and the product is
    converted back to the G-basis with ``to_groth_with_ring``.
    """

    def __init__(self, genset, coeff_genset, beta=None):
        super().__init__(genset, coeff_genset)
        if beta is None:
            beta = Symbol("\u03b2")
        self._beta = beta
        self._double_schubert_ring = DoubleSchubertRing(genset, coeff_genset)
        self.dtype = type("DoubleGrothendieckElement", (DoubleGrothendieckElement,), {"ring": self})

    def __str__(self):
        return f"Double Grothendieck polynomial ring in {self.genset.label} and {self.coeff_genset.label}"

    def __repr__(self):
        return f"DoubleGrothendieckRing({self.genset!r}, {self.coeff_genset!r}, beta={self._beta!r})"

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "DGR", self._beta))

    def __eq__(self, other):
        return type(self) is type(other) and self.genset == other.genset and self.coeff_genset == other.coeff_genset and self._beta == other._beta

    def perm_subs(self, elem, perm):
        elem_schub = self._as_schub(elem)
        dct = {self.genset[i]: -self.coeff_genset[perm[i - 1]]/(S.One + self._beta * self.coeff_genset[perm[i - 1]]) for i in range(1, max([len(p) for p in elem_schub.keys()]) + 1)}
        return elem_schub.eval(dct)

    @property
    def beta(self):
        return self._beta

    @cache
    def _as_schub_cached(self, perm):
        ring = self._double_schubert_ring
        return grothendieck_poly_with_ring(perm, ring, self._beta, keep_as_schub=True)

    def _as_schub(self, elem):
        """Expand a G-basis element into the underlying DoubleSchubertRing."""
        result = self._double_schubert_ring.zero
        for k, v in elem.items():
            result += v * self._as_schub_cached(k)
        return result

    @cached_property
    def vanish_subs_dict(self):
        ring = self._double_schubert_ring
        return {ring.genset[i]: -ring.coeff_genset[i] / (S.One + self._beta * ring.coeff_genset[i]) for i in range(50)}

    def permuted_subs_dict(self, perm):
        ring = self._double_schubert_ring
        return {ring.genset[i]: -ring.coeff_genset[perm[i - 1]] / (S.One + self._beta * ring.coeff_genset[perm[i - 1]]) for i in range(1, len(perm) + 1)}

    def from_double_schubert_elem(self, elem, _simplify=False):
        from schubmult.symbolic import expand

        ring = self._double_schubert_ring

        val = ring.from_dict({k: sympify(v).expand() for k, v in elem.items() if sympify(v).expand() != S.Zero})

        final_result = self.zero
        # beta = self._beta

        checked = set()
        last_val = val

        @cache
        def _isobaric_perm(perm1, div_perm, _beta):
            coeff = S.One
            if div_perm.inv == 0:
                return ring.from_dict({perm1: coeff})
            desc = min((~div_perm).descents())
            return _isobaric_perm(perm1, ~((~div_perm).swap(desc, desc + 1)), _beta).isobaric(desc + 1, _beta)

        def _isobaric_schub_elem(elem, div_perm, _beta):
            result = ring.zero
            for k, v in elem.items():
                result += v * _isobaric_perm(k, div_perm, _beta)
            return result


        while not val.almosteq(ring.zero):
            residual = S.Zero
            while True:
                val = ring.from_dict({k: v.expand() for k, v in val.items() if v.expand() != S.Zero})
                some_perm = min(set(val.keys()) - checked, key=lambda k: (k.inv, k), default=None)
                if some_perm is None:
                    if _simplify:
                        return self.from_dict({k: v.simplify() for k, v in final_result.items()})
                    return final_result
                residual = val.eval(self.permuted_subs_dict(some_perm)).simplify()
                checked.add(some_perm)

                if expand(residual) == S.Zero:
                    continue
                residual = (residual / self._as_schub_cached(some_perm).eval(self.permuted_subs_dict(some_perm))).simplify()
                break
            final_result += residual * self(some_perm)
            if some_perm.inv == 0:
                val = val - residual * ring.one
            else:
                val = val - residual * self._as_schub_cached(some_perm)
            if val.almosteq(last_val):
                raise ValueError(f"Failed to reduce {last_val} further; got stuck at {val}")
            last_val = val
        if _simplify:
            return self.from_dict({k: v.simplify() for k, v in final_result.items()})
        return final_result

    # def from_double_schubert_elem_fixed(self, elem, _simplify=False):
    #     from schubmult.symbolic import efficient_subs, expand, sympify
    #     from sympy import div

    #     ring = self._double_schubert_ring

    #     val = ring.from_dict({k: sympify(v).expand() for k, v in elem.items() if sympify(v).expand() != S.Zero})

    #     final_result = self.zero
    #     beta = self._beta

    #     fracto_subs = self.vanish_subs_dict
    #     checked = set()
    #     last_val = val

    #     while not val.almosteq(ring.zero):
    #         residual = S.Zero
    #         while True:
    #             val = ring.from_dict({k: v.expand() for k, v in val.items() if v.expand() != S.Zero})
    #             some_perm = min(set(val.keys()) - checked, key=lambda k: (k.inv, k), default=None)
    #             if some_perm is None:
    #                 if _simplify:
    #                     return self.from_dict({k: v.simplify() for k, v in final_result.items()})
    #                 return final_result
    #             if some_perm.inv == 0:
    #                 residual = efficient_subs(val.as_polynomial(), fracto_subs).simplify()
    #             else:
    #                 # iso_val = _isobaric_schub_elem(val, some_perm, beta)
    #                 fracto_wacto_subs = {ring.genset[(~some_perm)[i] - 1]: -ring.coeff_genset[i] / (S.One + beta * ring.coeff_genset[i]) for i in range(50)}
    #                 residual = efficient_subs(val.as_polynomial(), fracto_wacto_subs).simplify()
    #                 for (a,b) in some_perm.inversion_set:
    #                     residual, _ = div(residual, (ring.genset[b]))
    #             checked.add(some_perm)
    #             residual = expand(residual)
    #             if residual == S.Zero:
    #                 continue
    #             break
    #         final_result += residual * self(some_perm)
    #         if some_perm.inv == 0:
    #             val = val - residual * ring.one
    #         else:
    #             val = val - residual * self._as_schub_cached(some_perm)
    #         if val.almosteq(last_val):
    #             raise ValueError(f"Failed to reduce {last_val} further; got stuck at {val}")
    #         last_val = val
    #     if _simplify:
    #         return self.from_dict({k: v.simplify() for k, v in final_result.items()})
    #     return final_result

    def mul(self, elem, other):
        if not isinstance(other, BaseSchubertElement):
            return super().mul(elem, other)
        return self.from_double_schubert_elem(self._as_schub(elem) * other.ring._as_schub(other))

    def from_expr(self, expr):
        return self.from_double_schubert_elem(self._double_schubert_ring.from_expr(expr))

    @cache
    def cached_schubpoly(self, k):
        return grothendieck_poly_with_ring(k, self._double_schubert_ring, self._beta)

    def printing_term(self, k, prefix=""):
        return spolymod.DoubleGrothendieckPoly(k, self.genset.label, self.coeff_genset.label, prefix=prefix)

    def __call__(self, x):
        return self.new(x)

    def new(self, x):
        if isinstance(x, (list, tuple)):
            return self.from_dict({Permutation(x): S.One})
        if isinstance(x, Permutation):
            return self.from_dict({x: S.One})
        if isinstance(x, DoubleGrothendieckElement):
            if x.ring == self:
                return x
            raise ValueError("Different generating sets")
        return self.from_expr(x)

    def from_dict(self, dct):
        dct = {k: v for k, v in dct.items() if sympify(v).expand() != 0}
        return self.dtype(dct)


#DGx = DoubleGrothendieckRing(GeneratingSet("x"), GeneratingSet("y"))

def DGx(x, genset=GeneratingSet("y")):
    if isinstance(genset, str):
        genset = GeneratingSet(genset)
    return DoubleGrothendieckRing(GeneratingSet("x"), genset)(x)
