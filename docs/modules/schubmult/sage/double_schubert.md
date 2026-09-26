<a id="schubmult.sage.double_schubert"></a>

# schubmult.sage.double\_schubert

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

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomialRing"></a>

#### DoubleSchubertPolynomialRing

```python
def DoubleSchubertPolynomialRing(R,
                                 alphabet="y",
                                 coefficient_alphabets=("y", "z"))
```

Return the ring of double Schubert polynomials `\mathfrak{S}_w(x; \text{alphabet})` over ``R``.

INPUT:

- ``R`` -- a commutative ring (the scalars; the base ring of the result is the infinite
  polynomial ring ``R[alphabets]``)
- ``alphabet`` -- (default: ``'y'``) the letter of the second alphabet of the basis elements
- ``coefficient_alphabets`` -- (default: ``('y', 'z')``) letters available in coefficients;
  ``alphabet`` is always included. Rings over ``R`` with the same set of letters share a base
  ring, which is what lets elements of one coerce into another (mixed products).

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

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomial_class"></a>

## DoubleSchubertPolynomial\_class Objects

```python
class DoubleSchubertPolynomial_class(SchubmultBackedElement)
```

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomial_class.expand"></a>

#### expand

```python
def expand()
```

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

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomial_class.divided_difference"></a>

#### divided\_difference

```python
def divided_difference(i)
```

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

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomialRing_xbasis"></a>

## DoubleSchubertPolynomialRing\_xbasis Objects

```python
class DoubleSchubertPolynomialRing_xbasis(SchubmultBackedRing)
```

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomialRing_xbasis.__init__"></a>

#### \_\_init\_\_

```python
def __init__(R, alphabet, alphabets)
```

EXAMPLES::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: X = DoubleSchubertPolynomialRing(QQ)
    sage: X == loads(dumps(X))
    True
    sage: X is DoubleSchubertPolynomialRing(QQ, 'y', ('z', 'y'))
    True

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomialRing_xbasis.alphabet"></a>

#### alphabet

```python
def alphabet()
```

The letter of the second alphabet of the basis elements.

EXAMPLES::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: DoubleSchubertPolynomialRing(QQ, 'z').alphabet()
    'z'

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomialRing_xbasis.one_basis"></a>

#### one\_basis

```python
def one_basis()
```

EXAMPLES::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: DoubleSchubertPolynomialRing(QQ).one()
    X_y[1]

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomialRing_xbasis.degree_on_basis"></a>

#### degree\_on\_basis

```python
def degree_on_basis(w)
```

The degree of `\mathfrak{S}_w` is the length of `w`.

EXAMPLES::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: DoubleSchubertPolynomialRing(QQ)([3, 1, 2]).degree()
    2

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomialRing_xbasis.product_on_basis"></a>

#### product\_on\_basis

```python
def product_on_basis(left, right)
```

`\mathfrak{S}_u(x; y) \mathfrak{S}_v(x; y) = \sum_w c^w_{uv}(y) \mathfrak{S}_w(x; y)` via ``schubmult_double``.

EXAMPLES::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: X = DoubleSchubertPolynomialRing(QQ)
    sage: X.product_on_basis(Permutation([3, 2, 1]), Permutation([2, 1, 3]))
    (y_2-y_0)*X_y[3, 2, 1] + X_y[4, 2, 1, 3]

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomialRing_xbasis._element_constructor_"></a>

#### \_element\_constructor\_

```python
def _element_constructor_(x)
```

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

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomialRing_xbasis._coerce_map_from_"></a>

#### \_coerce\_map\_from\_

```python
def _coerce_map_from_(S)
```

Ordinary Schubert polynomial rings and double rings in another alphabet (over a base that
coerces into ours) coerce in.

EXAMPLES::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: X = DoubleSchubertPolynomialRing(QQ)
    sage: X.has_coerce_map_from(SchubertPolynomialRing(ZZ))
    True
    sage: X.has_coerce_map_from(DoubleSchubertPolynomialRing(QQ, 'z'))
    True
    sage: X.has_coerce_map_from(DoubleSchubertPolynomialRing(QQ, 'w'))
    False

<a id="schubmult.sage.double_schubert.DoubleSchubertPolynomialRing_xbasis.some_elements"></a>

#### some\_elements

```python
def some_elements()
```

EXAMPLES::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: DoubleSchubertPolynomialRing(QQ).some_elements()
    [X_y[1], X_y[1] + 2*X_y[2, 1], -X_y[3, 2, 1] + X_y[4, 2, 1, 3]]

