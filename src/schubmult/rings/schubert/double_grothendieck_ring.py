"""Double (equivariant K-theoretic) Grothendieck polynomial ring: the ``DGx`` interface.

`DoubleGrothendieckRing` represents ``G_w(x; y)`` with deformation parameter
``beta``. Products go through `schubmult.mult.groth_double.grothmult_double`
(the K-theoretic Monk/Pieri machinery) where available, with a fallback that
expands into the underlying `DoubleSchubertRing`. The ring also exposes the
localization/vanishing data used to convert between the Schubert and
Grothendieck bases (``permuted_subs_dict``, ``product_of_roots``,
``exp_root``, ``chevalley``).
"""

from functools import cache, cached_property

import schubmult.mult.groth as py
import schubmult.mult.groth_double as yz
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
        """Expand to an explicit polynomial: ``sum coeff * G_w(x; y)``."""
        from schubmult.symbolic import Add

        return Add(*[v * self.ring.cached_schubpoly(k) for k, v in self.items()])

    def perm_subs(self, perm):
        """Localize at the torus fixed point ``perm``: substitute ``x_i -> (-) y_{perm(i)}`` (formal inverse)."""
        return self.ring.perm_subs(self, perm)

    def simplify(self, factor=True):
        """Return a copy with each coefficient put in cancelled (and, by default, factored) rational
        normal form in ``y`` and ``beta``, dropping terms whose coefficient simplifies to zero.

        Products in this ring leave coefficients as unsimplified rational expressions; this
        makes them readable, e.g. ``(y_1 - y_2)/(1 + beta*y_2)``.
        """
        import sympy

        from schubmult.symbolic import sympify_sympy

        new_dict = {}
        for perm, coeff in self.items():
            expr = sympy.cancel(sympify_sympy(coeff))
            if expr == 0:
                continue
            if factor:
                expr = sympy.factor(expr)
            new_dict[perm] = sympify(expr)
        return self.ring.from_dict(new_dict)


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
        """Localize ``elem`` at ``perm``: expand into double Schubert polynomials and substitute
        ``x_i -> -y_{perm(i)} / (1 + beta y_{perm(i)})``.
        """
        elem_schub = self._as_schub(elem)
        dct = {self.genset[i]: -self.coeff_genset[perm[i - 1]] / (S.One + self._beta * self.coeff_genset[perm[i - 1]]) for i in range(1, max([len(p) for p in elem_schub.keys()]) + 1)}
        return elem_schub.eval(dct)

    @property
    def double_mul(self):
        """`schubmult.mult.groth_double.grothmult_double`."""
        return yz.grothmult_double

    @property
    def single_mul(self):
        """`schubmult.mult.groth.grothmult_py`."""
        return py.grothmult_py

    @property
    def beta(self):
        """The deformation parameter."""
        return self._beta

    @property
    def single_variable(self):
        """`schubmult.mult.groth_double.single_variable_groth` (K-theoretic Monk rule for ``x_k``)."""
        return yz.single_variable_groth

    @cache
    def _as_schub_cached(self, perm):
        """``G_perm`` expanded in the double Schubert basis (cached)."""
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
        """Localization at the identity: ``x_i -> -y_i / (1 + beta y_i)`` for the first 50 variables."""
        ring = self._double_schubert_ring
        return {ring.genset[i]: -ring.coeff_genset[i] / (S.One + self._beta * ring.coeff_genset[i]) for i in range(50)}

    def permuted_subs_dict(self, perm, length=None):
        """Localization at ``perm``: ``x_i -> -y_{perm(i)} / (1 + beta y_{perm(i)})`` for ``i <= length``."""
        ring = self._double_schubert_ring
        # perm[i - 1] returns i past the end of perm, so this stays correct for length > len(perm)
        if length is None:
            length = len(perm)
        return {ring.genset[i]: -ring.coeff_genset[perm[i - 1]] / (S.One + self._beta * ring.coeff_genset[perm[i - 1]]) for i in range(1, length + 1)}

    @cache
    def _weight(self, index):
        """The formal inverse ``(-) y_index = -y_index / (1 + beta y_index)``."""
        return -self.coeff_genset[index] / (S.One + self._beta * self.coeff_genset[index])

    @cache
    def product_of_roots(self, perm):
        """``prod_{(a,b) in Inv(perm^-1)} (y_a (-) y_b)``: the localization of ``G_perm`` at itself (Euler class)."""
        roots = (~perm).inversion_set
        result = S.One
        for a, b in roots:
            result *= (self.coeff_genset[a] - self.coeff_genset[b]) / (1 + self._beta * self.coeff_genset[b])
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

        from schubmult.symbolic import sympify_sympy

        return sympify(cancel(sympify_sympy(expr)))

    def div_by_product_of_roots(self, expr, perm):
        """Divide ``expr`` by ``product_of_roots(perm)`` one linear factor at a time, leaving any
        non-exact factors in the denominator (so the result may be a rational function).
        """
        from sympy import cancel, div, fraction


        num, den = fraction(cancel(expr))
        # num = expand(num)

        # Divide one linear factor at a time. Whatever does not divide exactly stays
        # in the denominator, so the result is a rational function in general.
        for a, b in sorted((~perm).inversion_set):
            quotient, remainder = div(num, (self.coeff_genset[a] - self.coeff_genset[b]))
            if remainder == S.Zero:
                num = quotient
            else:
                den = den * (self.coeff_genset[a] - self.coeff_genset[b])
            num = num * (1 + self._beta * self.coeff_genset[b])
        return sympify(cancel(num / den))

    @property
    def mult_poly_double(self):
        """`schubmult.mult.groth_double.mult_poly_groth_double`."""
        return yz.mult_poly_groth_double

    @property
    def mult_poly_single(self):
        """`schubmult.mult.groth.mult_poly_groth`."""
        return py.mult_poly_groth

    @cache
    def schub_as_groth(self, perm):
        """The double Schubert polynomial ``S_perm`` expanded in the Grothendieck basis (cached)."""
        return self._from_double_schubert_elem(self._double_schubert_ring(perm))

    def from_double_schubert_elem(self, elem):
        """Convert a `DoubleSchubertElement` into this ring's Grothendieck basis."""
        return self._from_double_schubert_elem(elem)
        # return sum(v * self.schub_as_groth(schub) for schub, v in elem.items())

    def _from_double_schubert_elem(self, elem):
        """Triangular basis change Schubert -> Grothendieck by repeated localization: pick the smallest
        remaining permutation, evaluate at its fixed point, divide by its Euler class to read off
        the ``G`` coefficient, subtract, and repeat.
        """
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
        """Multiply term by term via ``double_mul``; if a term raises ``NotImplementedError`` (the K-Monk
        rule doesn't cover it), fall back to multiplying by that term's explicit polynomial.
        """
        ring2 = elem2.ring
        result = self.zero

        for k, v in elem2.items():
            try:
                # raise NotImplementedError()
                result += v * self.from_dict(self.double_mul(elem, k, var2=self.coeff_genset, var3=ring2.coeff_genset, beta=self._beta))
            except NotImplementedError:
                # Fall back on the single basis element G_k, not on all of elem2.
                result += v * self.mul_expr(elem, ring2.cached_schubpoly(k))
        return result

    def mul_expr(self, elem, x):
        """Multiply by an expression: single ``x`` variables via the K-Monk rule, ``Add``/``Mul``/``Pow``
        recursively, anything else as a coefficient.
        """
        from schubmult.symbolic import Add, DomainElement, Mul, Pow

        if isinstance(x, DomainElement):
            raise TypeError(f"Cannot multiply {type(elem)} with {type(x)}")
        x = sympify(x)
        ind = self.genset.index(x)
        if ind != -1:
            return self.from_dict(self.single_variable(elem, ind, self.coeff_genset, beta=self._beta))
        if isinstance(x, Add):
            return self.sum([self.mul_expr(elem, arg) for arg in x.args])
        if isinstance(x, Mul):
            res = elem
            for arg in x.args:
                res = self.mul_expr(res, arg)
            return res
        if isinstance(x, Pow):
            res = elem
            base = x.args[0]
            exponent = x.args[1]
            if x.args[1] >= 0:
                for _ in range(int(exponent)):
                    res = self.mul_expr(res, base)
                return res
        return self.from_dict({k: v * self.domain_new(x) for k, v in elem.items()})

    def mul(self, elem, other):
        """Ring product via ``_best_effort_grothmult_double``."""
        # return self.from_double_schubert_elem(self._as_schub(elem) * other.ring._as_schub(other))
        return self._best_effort_grothmult_double(elem, other)

    def from_expr(self, expr):
        """Convert a polynomial into the Grothendieck basis by multiplying the identity by it."""
        return self.mul_expr(self.one, expr)
        # self.from_double_schubert_elem(self._double_schubert_ring.from_expr(expr))

    @cache
    def cached_schubpoly(self, k):
        """The explicit ``G_k(x; y)``, as a sum of `WCGraph` monomials weighted by ``beta^excess``."""
        from schubmult.combinatorics.wc_graph import WCGraph

        # return grothendieck_poly_with_ring(k, self._double_schubert_ring, self._beta)
        return sum([wc.polyvalue(self.genset, self.coeff_genset, beta=self._beta, prop_beta=True) for wc in WCGraph.all_wc_graphs(k)])

    def printing_term(self, k, prefix=""):
        """The ``DoubleGrothendieckPoly`` display symbol for basis element ``k``."""
        return spolymod.DoubleGrothendieckPoly(k, self.genset.label, self.coeff_genset.label, prefix=prefix)

    def __call__(self, x):
        return self.new(x)

    def new(self, x):
        """Build an element from a permutation/Lehmer list, an element of this ring, or a polynomial expression."""
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
        """Build an element from ``{Permutation: coeff}``, dropping exact zeros."""
        dct = {k: v for k, v in dct.items() if v != S.Zero}
        return self.dtype(dct)


# DGx = DoubleGrothendieckRing(GeneratingSet("x"), GeneratingSet("y"))


def DGx(x, genset=GeneratingSet("y")):
    """Construct a double Grothendieck element in ``x`` with coefficient alphabet ``genset`` (a
    `GeneratingSet`, a label string, or ``"0"`` for the zero alphabet).
    """
    from schubmult.symbolic.poly.variables import ZeroGeneratingSet

    if isinstance(genset, str):
        if genset == "0":
            genset = ZeroGeneratingSet()
        else:
            genset = GeneratingSet(genset)
    return DoubleGrothendieckRing(GeneratingSet("x"), genset)(x)
