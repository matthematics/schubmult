r"""
Double Schubert polynomials

The double Schubert polynomials `\mathfrak{S}_w(x; y)` are indexed by permutations `w` and form a
basis of `\ZZ[y_1, y_2, \ldots][x_1, x_2, \ldots]` over `\ZZ[y]`; they represent the equivariant
Schubert classes of the flag variety. Products are computed by the ``schubmult`` kernel
(:func:`schubmult.mult.double.schubmult_double`), so the structure constants come out as
polynomials in the second alphabet with no expansion to monomials.

Variables are 0-indexed on the Sage side, like :class:`~sage.combinat.schubert_polynomial.SchubertPolynomialRing`:
``expand()`` lands in ``x0, x1, ...`` and the coefficient alphabet is ``y_0, y_1, ...``. The sign
convention is `\mathfrak{S}_{21}(x; y) = x_0 - y_0`.

EXAMPLES::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: X = DoubleSchubertPolynomialRing(QQ); X
    Double Schubert polynomial ring in the alphabet y with X_y basis over Rational Field
    sage: X([3, 1, 2]) * X([2, 1])
    (y_2-y_0)*X_y[3, 1, 2] + X_y[4, 1, 2, 3]
    sage: X([3, 1, 2]).expand()
    x0^2 - x0*y0 - x0*y1 + y0*y1

Products agree with polynomial multiplication::

    sage: f = X([3, 1, 2]) * X([2, 1])
    sage: f.expand() == X([3, 1, 2]).expand() * X([2, 1]).expand()
    True

Setting the second alphabet to zero recovers ordinary Schubert polynomials::

    sage: S = SchubertPolynomialRing(QQ)
    sage: p = X([3, 1, 2]).expand()
    sage: S(p.subs({v: 0 for v in p.parent().gens()[3:]}))
    X[3, 1, 2]

Ordinary Schubert polynomials coerce in (they are polynomials in `x`, expanded in the double basis)::

    sage: X(S([2, 1]))
    y_0*X_y[1] + X_y[2, 1]
    sage: X([2, 1]) + S([2, 1])
    y_0*X_y[1] + 2*X_y[2, 1]

Mixed products `\mathfrak{S}_u(x; y) \mathfrak{S}_v(x; z)` expanded in the `y`-basis: build the
second factor in the ring with alphabet ``z``; it coerces into the ``y`` ring::

    sage: Z = DoubleSchubertPolynomialRing(QQ, 'z')
    sage: X([2, 1]) * Z([2, 1])
    (y_1-z_0)*X_y[2, 1] + X_y[3, 1, 2]

Arbitrary polynomials in `x` and the alphabets are expanded in the basis::

    sage: R.<x0, x1, y0, y1> = QQ[]
    sage: X(x0 - y0)
    X_y[2, 1]
    sage: X(x0*x1)
    y_1*y_0*X_y[1] + y_0*X_y[1, 3, 2] + X_y[2, 3, 1]
"""

from ._common import SchubmultBackedElement, SchubmultBackedRing, genset, to_sage_perm, to_schubmult_perm


def DoubleSchubertPolynomialRing(R, alphabet="y", coefficient_alphabets=("y", "z"), raw_coefficients=False):
    r"""
    Return the ring of double Schubert polynomials `\mathfrak{S}_w(x; \text{alphabet})` over ``R``.

    INPUT:

    - ``R`` -- a commutative ring (the scalars; the base ring of the result is the infinite
      polynomial ring ``R[alphabets]``)
    - ``alphabet`` -- (default: ``'y'``) the letter of the second alphabet of the basis elements
    - ``coefficient_alphabets`` -- (default: ``('y', 'z')``) letters available in coefficients;
      ``alphabet`` is always included. Rings over ``R`` with the same set of letters share a base
      ring, which is what lets elements of one coerce into another (mixed products).
    - ``raw_coefficients`` -- (default: ``False``) if ``True``, the base ring is the ring of unexpanded
      SymEngine expressions (:class:`~schubmult.sage.symengine_ring.SymEngineRing`) and the structure
      constants are kept in the factored form the kernel produces, with schubmult's 1-based variable
      names; ``R`` must be ``ZZ`` or ``QQ``

    EXAMPLES::

        sage: from schubmult.sage import DoubleSchubertPolynomialRing
        sage: X = DoubleSchubertPolynomialRing(ZZ); X
        Double Schubert polynomial ring in the alphabet y with X_y basis over Integer Ring
        sage: X.base_ring()
        Infinite polynomial ring in y, z over Integer Ring
        sage: TestSuite(X).run()
        sage: X(1)
        X_y[1]
        sage: X([1, 2, 3]) * X([2, 1, 3])
        X_y[2, 1]
        sage: X([2, 1, 3]) * X([2, 1, 3])
        (y_1-y_0)*X_y[2, 1] + X_y[3, 1, 2]
        sage: a = X([2, 1, 3]) + X([3, 1, 2, 4]); a^2
        (y_1-y_0)*X_y[2, 1] + (y_2^2-y_2*y_1-y_2*y_0+2*y_2+y_1*y_0-2*y_0+1)*X_y[3, 1, 2] + (y_3+y_2-y_1-y_0+2)*X_y[4, 1, 2, 3] + X_y[5, 1, 2, 3, 4]

    With unexpanded coefficients (note the 1-based variables)::

        sage: Y = DoubleSchubertPolynomialRing(ZZ, raw_coefficients=True); Y
        Double Schubert polynomial ring in the alphabet y with X_y basis over Integer Ring with unexpanded coefficients
        sage: TestSuite(Y).run()
        sage: f = Y([3, 1, 2]) * Y([3, 1, 2]) * Y([1, 3, 2]); f
        ((-y_1+y_3)*(-y_2+y_3)**2)*X_y[3, 1, 2] + ((-y_1+y_3)*(-y_2+y_3))*X_y[3, 2, 1] + ((-y_1+y_3)*(-y_2+y_3)+(-y_1-y_2+y_3+y_4)*(-y_2+y_4))*X_y[4, 1, 2, 3] + (-y_1-y_2+y_3+y_4)*X_y[4, 2, 1, 3] + (-y_1-2*y_2+y_3+y_4+y_5)*X_y[5, 1, 2, 3, 4] + X_y[5, 2, 1, 3, 4] + X_y[6, 1, 2, 3, 4, 5]
        sage: Y(X([3, 1, 2]) * X([3, 1, 2]) * X([1, 3, 2])) == f, X(f) == X([3, 1, 2]) * X([3, 1, 2]) * X([1, 3, 2])
        (True, True)
        sage: f.expand() == Y([3, 1, 2]).expand()^2 * Y([1, 3, 2]).expand()
        True
    """
    names = tuple(sorted({str(alphabet), *map(str, coefficient_alphabets)}))
    return DoubleSchubertPolynomialRing_xbasis(R, str(alphabet), names, raw_coefficients)


class DoubleSchubertPolynomial_class(SchubmultBackedElement):
    def expand(self):
        r"""
        Expand into a polynomial in ``x0, x1, ...`` and the coefficient alphabets ``y0, y1, ...``.

        EXAMPLES::

            sage: from schubmult.sage import DoubleSchubertPolynomialRing
            sage: X = DoubleSchubertPolynomialRing(ZZ)
            sage: X([2, 1]).expand()
            x0 - y0
            sage: X([1, 3, 2]).expand()
            x0 + x1 - y0 - y1
            sage: [X(p).expand() for p in Permutations(3)]
            [1, x0 + x1 - y0 - y1, x0 - y0, x0*x1 - x0*y0 - x1*y0 + y0^2, x0^2 - x0*y0 - x0*y1 + y0*y1, x0^2*x1 - x0^2*y0 - x0*x1*y0 + x0*y0^2 - x0*x1*y1 + x0*y0*y1 + x1*y0*y1 - y0^2*y1]
            sage: X([3, 1, 2]).expand().parent()
            Multivariate Polynomial Ring in x0, x1, x2, y0, y1 over Integer Ring
        """
        return super().expand()

    def divided_difference(self, i):
        r"""
        The divided difference `\partial_i` in the `x` variables: `\mathfrak{S}_w \mapsto \mathfrak{S}_{w s_i}`
        when `i` is a descent of `w`, and `0` otherwise; the coefficients are untouched.

        EXAMPLES::

            sage: from schubmult.sage import DoubleSchubertPolynomialRing
            sage: X = DoubleSchubertPolynomialRing(QQ)
            sage: X([3, 2, 1]).divided_difference(1)
            X_y[2, 3, 1]
            sage: X([3, 2, 1]).divided_difference(2)
            X_y[3, 1, 2]
            sage: X([3, 1, 2]).divided_difference(2)
            0
            sage: f = X([3, 1, 2]) * X([2, 1])
            sage: g = f.expand(); x0, x1 = g.parent().gens()[:2]
            sage: (g - g.subs({x0: x1, x1: x0})) // (x0 - x1) == f.divided_difference(1).expand()
            True
        """
        out = {}
        for w, c in self:
            if i in w.descents():
                lst = list(w)
                lst[i - 1], lst[i] = lst[i], lst[i - 1]
                out[to_sage_perm(lst)] = c
        return self.parent()._from_dict(out)


class DoubleSchubertPolynomialRing_xbasis(SchubmultBackedRing):
    Element = DoubleSchubertPolynomial_class

    def __init__(self, R, alphabet, alphabets, raw=False):
        """
        EXAMPLES::

            sage: from schubmult.sage import DoubleSchubertPolynomialRing
            sage: X = DoubleSchubertPolynomialRing(QQ)
            sage: X == loads(dumps(X))
            True
            sage: X is DoubleSchubertPolynomialRing(QQ, 'y', ('z', 'y'))
            True
        """
        self._alphabet = alphabet
        super().__init__(R, alphabets, prefix=f"X_{alphabet}", name=f"Double Schubert polynomial ring in the alphabet {alphabet} with X_{alphabet} basis", raw=raw)

    def alphabet(self):
        """
        The letter of the second alphabet of the basis elements.

        EXAMPLES::

            sage: from schubmult.sage import DoubleSchubertPolynomialRing
            sage: DoubleSchubertPolynomialRing(QQ, 'z').alphabet()
            'z'
        """
        return self._alphabet

    def _schub_ring(self, alphabet=None):
        from schubmult.rings.schubert.double_schubert_ring import DoubleSchubertRing

        return DoubleSchubertRing(genset("x"), genset(alphabet or self._alphabet))

    def _basis_polynomial(self, w):
        from schubmult.symbolic.poly.schub_poly import schubpoly

        return schubpoly(to_schubmult_perm(w), genset("x"), genset(self._alphabet))

    def one_basis(self):
        """
        EXAMPLES::

            sage: from schubmult.sage import DoubleSchubertPolynomialRing
            sage: DoubleSchubertPolynomialRing(QQ).one()
            X_y[1]
        """
        return super().one_basis()

    def degree_on_basis(self, w):
        r"""
        The degree of `\mathfrak{S}_w` is the length of `w`.

        EXAMPLES::

            sage: from schubmult.sage import DoubleSchubertPolynomialRing
            sage: DoubleSchubertPolynomialRing(QQ)([3, 1, 2]).degree()
            2
        """
        return super().degree_on_basis(w)

    def product_on_basis(self, left, right):
        r"""
        `\mathfrak{S}_u(x; y) \mathfrak{S}_v(x; y) = \sum_w c^w_{uv}(y) \mathfrak{S}_w(x; y)` via ``schubmult_double``.

        EXAMPLES::

            sage: from schubmult.sage import DoubleSchubertPolynomialRing
            sage: X = DoubleSchubertPolynomialRing(QQ)
            sage: X.product_on_basis(Permutation([3, 2, 1]), Permutation([2, 1, 3]))
            (y_2-y_0)*X_y[3, 2, 1] + X_y[4, 2, 1, 3]
        """
        from schubmult.mult.double import schubmult_double
        from schubmult.symbolic import S

        ys = genset(self._alphabet)
        return self._convert_dict(schubmult_double({to_schubmult_perm(left): S.One}, to_schubmult_perm(right), ys, ys))

    def _element_constructor_(self, x):
        """
        EXAMPLES::

            sage: from schubmult.sage import DoubleSchubertPolynomialRing
            sage: X = DoubleSchubertPolynomialRing(QQ)
            sage: X([2, 1, 3])
            X_y[2, 1]
            sage: X(Permutation([2, 1, 3]))
            X_y[2, 1]
            sage: X([])
            X_y[1]
            sage: X([1, 2, 1])
            Traceback (most recent call last):
            ...
            ValueError: the input [1, 2, 1] is not a valid permutation

            sage: R.<x0, x1, x2, y0, y1> = QQ[]
            sage: X(x0^2*x1)
            y_1*y_0^2*X_y[1] + y_0^2*X_y[1, 3, 2] + y_1*y_0*X_y[2, 1] + (y_1+y_0)*X_y[2, 3, 1] + y_0*X_y[3, 1, 2] + X_y[3, 2, 1]
            sage: X(X([3, 2, 1]).expand()) == X([3, 2, 1])
            True
            sage: S.<x> = InfinitePolynomialRing(QQ)
            sage: X(x[0]^2*x[1]) == X(x0^2*x1)
            True

        Ordinary Schubert polynomials and the other-alphabet double rings coerce::

            sage: X(SchubertPolynomialRing(QQ)([3, 1, 2]))
            y_0^2*X_y[1] + (y_1+y_0)*X_y[2, 1] + X_y[3, 1, 2]
            sage: Z = DoubleSchubertPolynomialRing(QQ, 'z')
            sage: X(Z([2, 1]))
            (y_0-z_0)*X_y[1] + X_y[2, 1]
            sage: Z(X(Z([3, 1, 2]))) == Z([3, 1, 2])
            True

        Key and atom polynomials (whose variables are read as ``x``) too; symmetric functions need the
        number of variables, see :meth:`from_symmetric_function`::

            sage: k = KeyPolynomials(QQ)
            sage: X(k([0, 2]))
            (y_1^2+y_1*y_0+y_0^2)*X_y[1] + (y_2+y_1+y_0)*X_y[1, 3, 2] + X_y[1, 4, 2, 3]
            sage: X([2, 1]) + k([1])
            y_0*X_y[1] + 2*X_y[2, 1]
            sage: X(SymmetricFunctions(QQ).s()[2, 1])
            Traceback (most recent call last):
            ...
            TypeError: a symmetric function needs a number of variables: use Double Schubert polynomial ring in the alphabet y with X_y basis over Rational Field.from_symmetric_function(f, n)
        """
        return super()._element_constructor_(x)

    def _coerce_map_from_(self, S):
        """
        Ordinary Schubert polynomial rings, key and atom polynomial rings, and double rings in another
        alphabet (over a base that coerces into ours) coerce in.

        EXAMPLES::

            sage: from schubmult.sage import DoubleSchubertPolynomialRing
            sage: X = DoubleSchubertPolynomialRing(QQ)
            sage: X.has_coerce_map_from(SchubertPolynomialRing(ZZ))
            True
            sage: X.has_coerce_map_from(KeyPolynomials(ZZ))
            True
            sage: X.has_coerce_map_from(DoubleSchubertPolynomialRing(QQ, 'z'))
            True
            sage: X.has_coerce_map_from(DoubleSchubertPolynomialRing(QQ, 'w'))
            False
        """
        return super()._coerce_map_from_(S)

    def some_elements(self):
        """
        EXAMPLES::

            sage: from schubmult.sage import DoubleSchubertPolynomialRing
            sage: DoubleSchubertPolynomialRing(QQ).some_elements()
            [X_y[1], X_y[1] + 2*X_y[2, 1], -X_y[3, 2, 1] + X_y[4, 2, 1, 3]]
        """
        return super().some_elements()
