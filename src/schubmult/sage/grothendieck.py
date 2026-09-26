r"""
Grothendieck polynomials

The `\beta`-Grothendieck polynomials `\mathfrak{G}^\beta_w(x)` (Fomin-Kirillov) represent Schubert
classes in the connective K-theory of the flag variety; `\beta = -1` gives the classical Grothendieck
polynomials of Lascoux-Schutzenberger and `\beta = 0` the Schubert polynomials. The double versions
`\mathfrak{G}^\beta_w(x; y)` are the equivariant classes; here the second alphabet enters through
`x \oplus y = x + y + \beta x y`, so `\mathfrak{G}_{21}(x; y) = x_0 + y_0 + \beta x_0 y_0`.

Products are computed by the schubmult kernels (:func:`schubmult.mult.groth.grothmult_py`,
:func:`schubmult.mult.groth_double.grothmult_double`). Structure constants of the single polynomials
are polynomials in `\beta`; those of the double polynomials are rational functions in `\beta` and the
second alphabet, so the double ring lives over the fraction field of ``R[beta][y, z]``.

Variables are 0-indexed on the Sage side (``x0, y_0``) as for the other rings in :mod:`schubmult.sage`;
the deformation parameter is ``beta``.

EXAMPLES::

    sage: from schubmult.sage import GrothendieckPolynomialRing, DoubleGrothendieckPolynomialRing
    sage: G = GrothendieckPolynomialRing(QQ); G
    Grothendieck polynomial ring with G basis over Rational Field
    sage: G([1, 3, 2]) * G([2, 1])
    G[2, 3, 1] + G[3, 1, 2] + beta*G[3, 2, 1]
    sage: G([1, 3, 2]).expand()
    x0*x1*beta + x0 + x1

At `\beta = 0` these are Schubert polynomials::

    sage: S = SchubertPolynomialRing(QQ)
    sage: f = G([1, 3, 2]) * G([2, 1]); p = f.expand()
    sage: S(p.subs({p.parent()('beta'): 0}))
    X[2, 3, 1] + X[3, 1, 2]

Double Grothendieck polynomials, with mixed alphabets as for the double Schubert ring::

    sage: GD = DoubleGrothendieckPolynomialRing(QQ); GD
    Double Grothendieck polynomial ring in the alphabet y with G_y basis over Rational Field
    sage: GD([2, 1]) * GD([2, 1])
    -((y_1-y_0)/(beta*y_1+1))*G_y[2, 1] + ((beta*y_0+1)/(beta*y_1+1))*G_y[3, 1, 2]
    sage: GD([2, 1]).expand()
    x0*y0*beta + x0 + y0
    sage: GZ = DoubleGrothendieckPolynomialRing(QQ, 'z')
    sage: GD([2, 1]) * GZ([2, 1])
    -((y_1-z_0)/(beta*y_1+1))*G_y[2, 1] + ((beta*z_0+1)/(beta*y_1+1))*G_y[3, 1, 2]

Products agree with polynomial multiplication (in the fraction field, for the double ring)::

    sage: f = GD([3, 1, 2]) * GD([2, 3, 1])
    sage: f.expand() == GD([3, 1, 2]).expand() * GD([2, 3, 1]).expand()
    True
"""

from sage.rings.fraction_field import FractionField_generic
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.richcmp import op_EQ, op_NE
from sage.structure.unique_representation import UniqueRepresentation

from ._common import SchubmultBackedElement, SchubmultBackedRing, genset

BETA = "\u03b2"  # schubmult's symbol for the deformation parameter
BETA_VARIABLE = "beta_0"  # its name inside the infinite polynomial ring (which only has indexed variables)


def GrothendieckPolynomialRing(R):
    r"""
    Return the ring of `\beta`-Grothendieck polynomials `\mathfrak{G}^\beta_w(x)` over ``R``.

    The base ring is ``R[beta]``. Ordinary Schubert polynomials (elements of Sage's
    :func:`SchubertPolynomialRing`) and polynomials in ``x0, x1, ...`` (and ``beta``) coerce in.

    EXAMPLES::

        sage: from schubmult.sage import GrothendieckPolynomialRing
        sage: G = GrothendieckPolynomialRing(ZZ); G
        Grothendieck polynomial ring with G basis over Integer Ring
        sage: G.base_ring()
        Univariate Polynomial Ring in beta over Integer Ring
        sage: TestSuite(G).run()
        sage: G([2, 1]) * G([2, 1])
        G[3, 1, 2]
        sage: G([2, 1]) * G([1, 3, 2])
        G[2, 3, 1] + G[3, 1, 2] + beta*G[3, 2, 1]
        sage: G(SchubertPolynomialRing(ZZ)([1, 3, 2]))
        G[1, 3, 2] - beta*G[2, 3, 1]
        sage: R.<x0, x1, beta> = ZZ[]
        sage: G(x0 + x1 + beta*x0*x1)
        G[1, 3, 2]
    """
    return GrothendieckPolynomialRing_gbasis(R)


def DoubleGrothendieckPolynomialRing(R, alphabet="y", coefficient_alphabets=("y", "z")):
    r"""
    Return the ring of double `\beta`-Grothendieck polynomials `\mathfrak{G}^\beta_w(x; \text{alphabet})` over ``R``.

    The base ring is the fraction field of ``R[beta][alphabets]`` (see :class:`GrothendieckCoefficientField`):
    equivariant K-theoretic structure constants are rational in `\beta` and the second alphabet, with
    denominators that are products of `1 + \beta y_i`. As for
    :func:`~schubmult.sage.DoubleSchubertPolynomialRing`, rings over ``R`` with the same set of letters
    share a base ring and coerce into each other (mixed products).

    EXAMPLES::

        sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
        sage: GD = DoubleGrothendieckPolynomialRing(QQ); GD
        Double Grothendieck polynomial ring in the alphabet y with G_y basis over Rational Field
        sage: GD.base_ring()
        Fraction Field of Infinite polynomial ring in beta, y, z over Rational Field
        sage: TestSuite(GD).run()
        sage: GD([2, 1]) * GD([1, 3, 2])
        G_y[2, 3, 1] + G_y[3, 1, 2] + beta*G_y[3, 2, 1]
        sage: GD([3, 1, 2]) * GD([2, 1])
        -((y_2-y_0)/(beta*y_2+1))*G_y[3, 1, 2] + ((beta*y_0+1)/(beta*y_2+1))*G_y[4, 1, 2, 3]
        sage: GD([2, 1]).expand()
        x0*y0*beta + x0 + y0

    Setting `\beta = 0` gives the double Schubert polynomial in `x` and `-y` (the second alphabet enters
    through `x \oplus y`, whereas :func:`~schubmult.sage.DoubleSchubertPolynomialRing` uses `x - y`)::

        sage: from schubmult.sage import DoubleSchubertPolynomialRing
        sage: p = GD([3, 1, 2]).expand(); T = p.parent()
        sage: q = DoubleSchubertPolynomialRing(QQ)([3, 1, 2]).expand()
        sage: p.subs({T('beta'): 0}) == T(q.subs({g: -g for g in q.parent().gens() if str(g).startswith('y')}))
        True
    """
    names = tuple(sorted({str(alphabet), *map(str, coefficient_alphabets)}))
    return DoubleGrothendieckPolynomialRing_gbasis(R, str(alphabet), names)


class GrothendieckPolynomial_class(SchubmultBackedElement):
    pass


class GrothendieckCoefficient(FractionFieldElement):
    """
    Element of :class:`GrothendieckCoefficientField`: the deformation parameter prints as ``beta``.

    EXAMPLES::

        sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
        sage: GD = DoubleGrothendieckPolynomialRing(QQ); b = GD.beta(); y = GD.base_ring().ring().gen(1)
        sage: (y[0] - y[1]) / (1 + b*y[1])
        (-y_1 + y_0)/(beta*y_1 + 1)
        sage: latex(_)
        \frac{-y_{1} + y_{0}}{\beta y_{1} + 1}
    """

    def _repr_(self):
        return super()._repr_().replace(BETA_VARIABLE, "beta")

    def _latex_(self):
        return super()._latex_().replace(r"\beta_{0}", r"\beta")

    def _richcmp_(self, other, op):
        # a/d == b/d iff a == b: skip the cross-multiplication (two products of big polynomials) that
        # the generic comparison does when the denominators agree, which they do for equal
        # structure constants computed twice (associativity checks, comparing products)
        if op in (op_EQ, op_NE) and self.denominator() == other.denominator():
            eq = self.numerator() == other.numerator()
            return eq if op == op_EQ else not eq
        return super()._richcmp_(other, op)


class GrothendieckCoefficientField(UniqueRepresentation, FractionField_generic):
    r"""
    The coefficient field ``Frac(R[beta, y_0, y_1, ..., z_0, ...])`` of a double Grothendieck ring.

    It is the fraction field of an ``InfinitePolynomialRing`` in which ``beta`` is the indexed
    variable ``beta_0`` -- that keeps every finite ring underneath a libsingular ring over ``R``
    (``R[beta][y_0, ...]`` would be the generic, very slow implementation, and gcds over ``R(beta)``
    are slow too), which is what makes coefficient arithmetic fast. Elements print ``beta_0`` as ``beta``.

    EXAMPLES::

        sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
        sage: F = DoubleGrothendieckPolynomialRing(QQ).base_ring(); F
        Fraction Field of Infinite polynomial ring in beta, y, z over Rational Field
        sage: F is loads(dumps(F))
        True
        sage: F.ring().gens()
        (beta_*, y_*, z_*)
    """

    Element = GrothendieckCoefficient

    def __init__(self, R):
        FractionField_generic.__init__(self, R)
        self._element_class = self.element_class  # ``Element`` glued to the category's element class


class _GrothendieckMixin:
    _named_symbols = {BETA: "beta"}  # noqa: RUF012

    def beta(self):
        r"""
        The deformation parameter `\beta` as an element of the base ring.

        EXAMPLES::

            sage: from schubmult.sage import GrothendieckPolynomialRing
            sage: G = GrothendieckPolynomialRing(QQ)
            sage: G.beta() * G([2, 1])
            beta*G[2, 1]
        """
        return self._beta

    def _named_base_elements(self):
        return {BETA: self._beta_scalar}  # in a ring every underlying finite ring converts from

    def _schub_beta(self):
        from schubmult.symbolic import Symbol

        return Symbol(BETA)


class GrothendieckPolynomialRing_gbasis(_GrothendieckMixin, SchubmultBackedRing):
    Element = GrothendieckPolynomial_class

    def __init__(self, R):
        """
        EXAMPLES::

            sage: from schubmult.sage import GrothendieckPolynomialRing
            sage: G = GrothendieckPolynomialRing(QQ)
            sage: G == loads(dumps(G))
            True
        """
        self._alphabet = None
        super().__init__(R, (), prefix="G", name="Grothendieck polynomial ring with G basis")
        self._beta_scalar = self.base_ring().gen()
        self._beta = self._beta_scalar

    @staticmethod
    def _make_base_ring(R, alphabets):  # noqa: ARG004
        return PolynomialRing(R, "beta")

    def _schub_ring(self, alphabet=None):  # noqa: ARG002
        from schubmult.rings.schubert.grothendieck_ring import GrothendieckRing

        return GrothendieckRing(genset("x"), self._schub_beta())

    def product_on_basis(self, left, right):
        r"""
        `\mathfrak{G}_u \mathfrak{G}_v = \sum_w c^w_{uv}(\beta) \mathfrak{G}_w` via ``grothmult_py``.

        EXAMPLES::

            sage: from schubmult.sage import GrothendieckPolynomialRing
            sage: G = GrothendieckPolynomialRing(QQ)
            sage: G.product_on_basis(Permutation([1, 3, 2]), Permutation([1, 3, 2]))
            G[1, 4, 2, 3] + G[2, 3, 1] + beta*G[2, 4, 1, 3]
        """
        return super().product_on_basis(left, right)


class DoubleGrothendieckPolynomialRing_gbasis(_GrothendieckMixin, SchubmultBackedRing):
    Element = GrothendieckPolynomial_class
    _base_aliases = {BETA_VARIABLE: "beta"}  # noqa: RUF012

    def __init__(self, R, alphabet, alphabets):
        """
        EXAMPLES::

            sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
            sage: GD = DoubleGrothendieckPolynomialRing(QQ, 'z')
            sage: GD == loads(dumps(GD))
            True
            sage: GD is DoubleGrothendieckPolynomialRing(QQ, 'z', ('y',))
            True
        """
        self._alphabet = alphabet
        super().__init__(R, alphabets, prefix=f"G_{alphabet}", name=f"Double Grothendieck polynomial ring in the alphabet {alphabet} with G_{alphabet} basis")
        D = self.base_ring().ring()
        self._beta_scalar = D.gen(D.variable_names().index("beta"))[0].polynomial()
        self._beta = self.base_ring()(self._beta_scalar)

    @staticmethod
    def _make_base_ring(R, alphabets):
        return GrothendieckCoefficientField(InfinitePolynomialRing(R, ["beta", *alphabets]))

    def alphabet(self):
        """
        The letter of the second alphabet of the basis elements.

        EXAMPLES::

            sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
            sage: DoubleGrothendieckPolynomialRing(QQ, 'z').alphabet()
            'z'
        """
        return self._alphabet

    def _schub_ring(self, alphabet=None):
        from schubmult.rings.schubert.double_grothendieck_ring import DoubleGrothendieckRing

        return DoubleGrothendieckRing(genset("x"), genset(alphabet or self._alphabet), self._schub_beta())

    def some_elements(self):
        """
        A few small elements (structure constants grow quickly with the permutations, so these stay in `S_3`).

        EXAMPLES::

            sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
            sage: DoubleGrothendieckPolynomialRing(QQ).some_elements()
            [G_y[1], G_y[1] + 2*G_y[2, 1], -G_y[1, 3, 2] + G_y[2, 3, 1]]
        """
        return [self.one(), self.one() + 2 * self([2, 1]), self([2, 3, 1]) - self([1, 3, 2])]

    def product_on_basis(self, left, right):
        r"""
        `\mathfrak{G}_u(x; y) \mathfrak{G}_v(x; y) = \sum_w c^w_{uv}(\beta; y) \mathfrak{G}_w(x; y)` via ``grothmult_double``.

        EXAMPLES::

            sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
            sage: GD = DoubleGrothendieckPolynomialRing(QQ)
            sage: GD.product_on_basis(Permutation([2, 1]), Permutation([2, 1]))
            -((y_1-y_0)/(beta*y_1+1))*G_y[2, 1] + ((beta*y_0+1)/(beta*y_1+1))*G_y[3, 1, 2]
        """
        return super().product_on_basis(left, right)
