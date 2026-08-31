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

    def permuted_subs_dict(self, perm, length=None):
        ring = self._double_schubert_ring
        # perm[i - 1] returns i past the end of perm, so this stays correct for length > len(perm)
        if length is None:
            length = len(perm)
        return {ring.genset[i]: -ring.coeff_genset[perm[i - 1]] / (S.One + self._beta * ring.coeff_genset[perm[i - 1]]) for i in range(1, length + 1)}

    @cache
    def _weight(self, index):
        return -self.coeff_genset[index]/(S.One + self._beta * self.coeff_genset[index])

    @cache
    def product_of_roots(self, perm):
        roots = (~perm).inversion_set
        result = S.One
        for (a, b) in roots:
            result *= (self.coeff_genset[a] - self.coeff_genset[b])/(1+self._beta*self.coeff_genset[b])
        return result

    @cache
    def exp_root(self, a, b):
        """Polynomial form of the root ``eps_a - eps_b``, i.e. ``(x^root - 1)/beta``."""
        return (self.coeff_genset[a] - self.coeff_genset[b]) / (S.One + self._beta * self.coeff_genset[b])

    @cache
    def exp_weight(self, index):
        """Polynomial form of the fundamental weight ``eps_index``."""
        return self.coeff_genset[index] / (S.One + self._beta * self.coeff_genset[index])

    def exp_character(self, weight):
        """Polynomial form of the multiplicative character ``x^weight``.

        Determined by ``exp_root``: ``x^{eps_i} = 1 + beta * y_i``, so that
        ``(x^{eps_a - eps_b} - 1)/beta`` reproduces ``exp_root(a, b)``.
        """
        result = S.One
        for index, power in enumerate(weight, start=1):
            if power:
                result *= (S.One + self._beta * self.coeff_genset[index]) ** power
        return result

    def chevalley(self, weight, perm, n=None):
        """Lenart--Postnikov K_T-Chevalley formula for ``e^weight * G_perm``.

        ``weight`` is an integer vector in the ``epsilon`` basis.
        """
        from .chevalley import kt_chevalley_coefficients

        result = self.zero
        for (w, mu), coeff in kt_chevalley_coefficients(perm, weight, n).items():
            result += coeff * self.exp_character(mu) * self(w)
        return result

    def _cancel(self, expr):
        """Reduce a rational function to lowest terms without the overhead of full simplify()."""
        from sympy import cancel

        from schubmult.symbolic import sympify, sympify_sympy
        return sympify(cancel(sympify_sympy(expr)))

    def div_by_product_of_roots(self, expr, perm):
        from sympy import cancel, div, fraction

        from schubmult.symbolic import sympify

        num, den = fraction(cancel(expr))
        #num = expand(num)

        # Divide one linear factor at a time. Whatever does not divide exactly stays
        # in the denominator, so the result is a rational function in general.
        for a, b in sorted((~perm).inversion_set):
            quotient, remainder = div(num, (self.coeff_genset[a] - self.coeff_genset[b]))
            if remainder == S.Zero:
                num = quotient
            else:
                den = den * (self.coeff_genset[a] - self.coeff_genset[b])
            num = num * (1 + self._beta * self.coeff_genset[b])
        return sympify(cancel(num/den))


    @cache
    def schub_as_groth(self, perm):
        return self._from_double_schubert_elem(self._double_schubert_ring(perm))

    def from_double_schubert_elem(self, elem):
        return self._from_double_schubert_elem(elem)
        #return sum(v * self.schub_as_groth(schub) for schub, v in elem.items())

    def _from_double_schubert_elem(self, elem):
        from sympy import cancel

        from schubmult.symbolic import expand

        from .double_schubert_ring import DoubleSchubertElement

        ring = self._double_schubert_ring

        val = ring.from_dict({k: expand(cancel(v)) for k, v in elem.items()})
        val = val.strip_zeros()

        final_result = self.zero
        # beta = self._beta

        checked = set()
        # last_val = val

        while not val.almosteq(ring.zero):
            residual = S.Zero
            some_perm = min(set(val.keys()) - checked, key=lambda k: (k.inv, k), default=None)
            if some_perm is None:
                return final_result
            num_vars = max([*[len(k) for k in val.keys()], 0])
            residual = val.eval(self.permuted_subs_dict(some_perm, num_vars))
            if isinstance(residual, DoubleSchubertElement):
                residual = residual.as_polynomial()
            # residual is a rational function of coeff-genset/beta; cancel() collapses
            # it to lowest terms far faster than a full simplify() (which additionally
            # tries powsimp/trigsimp/etc. that never help here).
            # residual = self._cancel(residual)
            checked.add(some_perm)
            if expand(cancel(residual)) == S.Zero:
                continue
            residual = self.div_by_product_of_roots(residual, some_perm)
            final_result += residual * self(some_perm)
            if some_perm.inv == 0:
                val = val - residual * ring.one
            else:
                val = val - residual * self._as_schub_cached(some_perm)
            val = ring.from_dict({k: expand(cancel(v)) for k, v in val.items()})
            val = val.strip_zeros()
        return final_result

    def _best_effort_grothmult_double(self, elem, elem2):
        from schubmult.mult.groth_double import grothmult_double
        ring2 = elem2.ring
        result = self.zero

        for k, v in elem2.items():
            try:
                result += v * self.from_dict(grothmult_double(elem, k, var2=self.coeff_genset, var3=ring2.coeff_genset, beta=self._beta))
            except NotImplementedError:
                result += v * self.from_double_schubert_elem(self._as_schub(elem) * ring2._as_schub(elem2))
        return result


    def mul(self, elem, other):
        #return self.from_double_schubert_elem(self._as_schub(elem) * other.ring._as_schub(other))
        return self._best_effort_grothmult_double(elem, other)

    def from_expr(self, expr):
        return self.from_double_schubert_elem(self._double_schubert_ring.from_expr(expr))

    def mul_expr(self, elem, expr):
        return self.mul(elem, self.from_expr(expr))

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
    from schubmult.symbolic.poly.variables import ZeroGeneratingSet
    if isinstance(genset, str):
        if genset == "0":
            genset = ZeroGeneratingSet()
        else:
            genset = GeneratingSet(genset)
    return DoubleGrothendieckRing(GeneratingSet("x"), genset)(x)
