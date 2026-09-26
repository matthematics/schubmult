<a id="schubmult.sage.grothendieck"></a>

# schubmult.sage.grothendieck

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

<a id="schubmult.sage.grothendieck.BETA"></a>

#### BETA

schubmult's symbol for the deformation parameter

<a id="schubmult.sage.grothendieck.BETA_VARIABLE"></a>

#### BETA\_VARIABLE

its name inside the infinite polynomial ring (which only has indexed variables)

<a id="schubmult.sage.grothendieck.GrothendieckPolynomialRing"></a>

#### GrothendieckPolynomialRing

```python
def GrothendieckPolynomialRing(R)
```

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

<a id="schubmult.sage.grothendieck.DoubleGrothendieckPolynomialRing"></a>

#### DoubleGrothendieckPolynomialRing

```python
def DoubleGrothendieckPolynomialRing(R,
                                     alphabet="y",
                                     coefficient_alphabets=("y", "z"))
```

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

<a id="schubmult.sage.grothendieck.GrothendieckCoefficient"></a>

## GrothendieckCoefficient Objects

```python
class GrothendieckCoefficient(FractionFieldElement)
```

Element of :class:`GrothendieckCoefficientField`: the deformation parameter prints as ``beta``.

EXAMPLES::

    sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
    sage: GD = DoubleGrothendieckPolynomialRing(QQ); b = GD.beta(); y = GD.base_ring().ring().gen(1)
    sage: (y[0] - y[1]) / (1 + b*y[1])
    (-y_1 + y_0)/(beta*y_1 + 1)
    sage: latex(_)
    rac{-y_{1} + y_{0}}{eta y_{1} + 1}

<a id="schubmult.sage.grothendieck.GrothendieckCoefficientField"></a>

## GrothendieckCoefficientField Objects

```python
class GrothendieckCoefficientField(UniqueRepresentation,
                                   FractionField_generic)
```

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

<a id="schubmult.sage.grothendieck._GrothendieckMixin"></a>

## \_GrothendieckMixin Objects

```python
class _GrothendieckMixin()
```

<a id="schubmult.sage.grothendieck._GrothendieckMixin.beta"></a>

#### beta

```python
def beta()
```

The deformation parameter `\beta` as an element of the base ring.

EXAMPLES::

    sage: from schubmult.sage import GrothendieckPolynomialRing
    sage: G = GrothendieckPolynomialRing(QQ)
    sage: G.beta() * G([2, 1])
    beta*G[2, 1]

<a id="schubmult.sage.grothendieck.GrothendieckPolynomialRing_gbasis"></a>

## GrothendieckPolynomialRing\_gbasis Objects

```python
class GrothendieckPolynomialRing_gbasis(_GrothendieckMixin,
                                        SchubmultBackedRing)
```

<a id="schubmult.sage.grothendieck.GrothendieckPolynomialRing_gbasis.__init__"></a>

#### \_\_init\_\_

```python
def __init__(R)
```

EXAMPLES::

    sage: from schubmult.sage import GrothendieckPolynomialRing
    sage: G = GrothendieckPolynomialRing(QQ)
    sage: G == loads(dumps(G))
    True

<a id="schubmult.sage.grothendieck.GrothendieckPolynomialRing_gbasis.product_on_basis"></a>

#### product\_on\_basis

```python
def product_on_basis(left, right)
```

`\mathfrak{G}_u \mathfrak{G}_v = \sum_w c^w_{uv}(\beta) \mathfrak{G}_w` via ``grothmult_py``.

EXAMPLES::

    sage: from schubmult.sage import GrothendieckPolynomialRing
    sage: G = GrothendieckPolynomialRing(QQ)
    sage: G.product_on_basis(Permutation([1, 3, 2]), Permutation([1, 3, 2]))
    G[1, 4, 2, 3] + G[2, 3, 1] + beta*G[2, 4, 1, 3]

<a id="schubmult.sage.grothendieck.DoubleGrothendieckPolynomialRing_gbasis"></a>

## DoubleGrothendieckPolynomialRing\_gbasis Objects

```python
class DoubleGrothendieckPolynomialRing_gbasis(_GrothendieckMixin,
                                              SchubmultBackedRing)
```

<a id="schubmult.sage.grothendieck.DoubleGrothendieckPolynomialRing_gbasis.__init__"></a>

#### \_\_init\_\_

```python
def __init__(R, alphabet, alphabets)
```

EXAMPLES::

    sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
    sage: GD = DoubleGrothendieckPolynomialRing(QQ, 'z')
    sage: GD == loads(dumps(GD))
    True
    sage: GD is DoubleGrothendieckPolynomialRing(QQ, 'z', ('y',))
    True

<a id="schubmult.sage.grothendieck.DoubleGrothendieckPolynomialRing_gbasis.alphabet"></a>

#### alphabet

```python
def alphabet()
```

The letter of the second alphabet of the basis elements.

EXAMPLES::

    sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
    sage: DoubleGrothendieckPolynomialRing(QQ, 'z').alphabet()
    'z'

<a id="schubmult.sage.grothendieck.DoubleGrothendieckPolynomialRing_gbasis.some_elements"></a>

#### some\_elements

```python
def some_elements()
```

A few small elements (structure constants grow quickly with the permutations, so these stay in `S_3`).

EXAMPLES::

    sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
    sage: DoubleGrothendieckPolynomialRing(QQ).some_elements()
    [G_y[1], G_y[1] + 2*G_y[2, 1], -G_y[1, 3, 2] + G_y[2, 3, 1]]

<a id="schubmult.sage.grothendieck.DoubleGrothendieckPolynomialRing_gbasis.product_on_basis"></a>

#### product\_on\_basis

```python
def product_on_basis(left, right)
```

`\mathfrak{G}_u(x; y) \mathfrak{G}_v(x; y) = \sum_w c^w_{uv}(\beta; y) \mathfrak{G}_w(x; y)` via ``grothmult_double``.

EXAMPLES::

    sage: from schubmult.sage import DoubleGrothendieckPolynomialRing
    sage: GD = DoubleGrothendieckPolynomialRing(QQ)
    sage: GD.product_on_basis(Permutation([2, 1]), Permutation([2, 1]))
    -((y_1-y_0)/(beta*y_1+1))*G_y[2, 1] + ((beta*y_0+1)/(beta*y_1+1))*G_y[3, 1, 2]

