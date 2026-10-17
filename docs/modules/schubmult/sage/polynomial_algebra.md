<a id="schubmult.sage.polynomial_algebra"></a>

# schubmult.sage.polynomial\_algebra

The polynomial ring with its combinatorial bases

:func:`PolynomialAlgebra` is the polynomial ring `R[x_0, x_1, \ldots]` (or `R[x_0, \ldots, x_{n-1}]`)
as a parent with realizations (in the sense of
:class:`sage.categories.with_realizations.WithRealizations`, like :class:`SymmetricFunctions`): one
:class:`CombinatorialFreeModule` per basis, with coercions between them. The bases and their
conversions are computed by schubmult's :class:`~schubmult.rings.polynomial_algebra.PolynomialAlgebra`
bases and Schubert kernels.

Each basis is indexed the way its polynomials are indexed in the literature:

- ``schubert``, ``grothendieck``: by permutations, `\mathfrak S_w` and the Grothendieck polynomials
  `\mathfrak G_w` (prefixes ``S``, ``G``);
- ``monomial``, ``key``, ``fundamental_slide``, ``monomial_slide``, ``forest``, ``glide``, ``lascoux``,
  ``grove``: by weak compositions `\alpha`: `x^\alpha`, the key polynomials `\kappa_\alpha`, the slide
  polynomials of Assaf-Searles, the forest polynomials of Nadeau-Spink-Tewari, and the K-theoretic
  glide, Lascoux and grove polynomials (prefixes ``x``, ``k``, ``F``, ``M``, ``P``, ``Gl``, ``L``, ``Gr``);
- ``elementary``: products `\prod_{j < n} e_{a_j}(x_0, \ldots, x_{j-1}) \cdot \prod_i e_{b_i}(x_0, \ldots, x_{n-1})`
  of elementary symmetric polynomials, a basis of the polynomial ring in `n` variables *for each fixed
  `n`* -- so it is only available in ``PolynomialAlgebra(R, n)`` (prefix ``E``; see :meth:`PolynomialAlgebra.elementary`).

The K-theoretic bases (``grothendieck``, ``glide``, ``lascoux``, ``grove``) are taken at `\beta = -1`, the
classical convention (`\mathfrak G_{132} = x_0 + x_1 - x_0 x_1`). Nothing is lost: with `\deg \beta = -1`
these polynomials are homogeneous, so `P^{\beta}_a(x) = (-\beta)^{|a|} P^{-1}_a(-x/\beta)` recovers any
`\beta`; schubmult's own bases work at `\beta = 1` and the layer applies that sign twist.

schubmult's ``PolynomialAlgebra`` is graded by the number of variables (a key of length `n` lives in
the `n`-variable slice and keys of different lengths multiply to zero -- the structure dual to the
free algebra). ``PolynomialAlgebra(R)`` is the plain polynomial ring instead: trailing zeros of a
composition do not matter, `\mathfrak S_w` is `\mathfrak S_w` however many variables are in play, and
products are polynomial products; the layer pads keys to a common number of variables before calling
schubmult. ``PolynomialAlgebra(R, n)`` is the `n`-variable slice itself, where the number of variables
is part of the indexing (compositions have length `n`, permutations have their last descent at most
`n`, and the elementary basis is available).

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: A = PolynomialAlgebra(QQ); A
    Polynomial ring in x0, x1, ... over Rational Field with combinatorial bases
    sage: S = A.schubert(); k = A.key(); x = A.monomial()
    sage: S[3, 1, 2]
    S[3, 1, 2]
    sage: k(S[3, 1, 2])
    k[2]
    sage: x(k[2, 0, 1])
    x[2, 0, 1] + x[2, 1]
    sage: k[2, 0, 1].expand()
    x0^2*x1 + x0^2*x2
    sage: S[2, 1, 4, 3] * S[1, 3, 2]
    S[2, 3, 4, 1] + S[2, 4, 1, 3] + S[3, 1, 4, 2] + S[4, 1, 2, 3]
    sage: A.grothendieck()(S[1, 3, 2])
    G[1, 3, 2] + G[2, 3, 1]

The bases coerce into one another, and Sage's own Schubert and key polynomials coerce in::

    sage: k[1, 0, 2] + S[1, 3, 2]
    k[0, 1] + k[1, 0, 2]
    sage: A.forest()(SchubertPolynomialRing(QQ)([3, 1, 2]))
    P[2]
    sage: A.fundamental_slide()(KeyPolynomials(QQ)([1, 0, 2]))
    F[1, 0, 2] + F[2, 0, 1]

In a fixed number of variables::

    sage: B = PolynomialAlgebra(QQ, 3); B
    Polynomial ring in x0, x1, x2 over Rational Field with combinatorial bases
    sage: E = B.elementary(); E[1, 0, 2]
    E[1, 0, 2]
    sage: E[1, 0, 2].expand()
    x0^2*x1 + x0^2*x2 + x0*x1*x2
    sage: B.key()(E[1, 0, 2])
    k[1, 1, 1] + k[2, 0, 1]
    sage: B.key()[2, 0, 1] * B.key()[0, 1, 0]
    k[2, 1, 1] + k[2, 2, 0] + k[3, 0, 1]

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra"></a>

## PolynomialAlgebra Objects

```python
class PolynomialAlgebra(UniqueRepresentation, Parent)
```

The polynomial ring `R[x_0, x_1, \ldots]`, or `R[x_0, \ldots, x_{n-1}]` if ``n`` is given, with its
combinatorial bases as realizations.

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: A = PolynomialAlgebra(ZZ); A
    Polynomial ring in x0, x1, ... over Integer Ring with combinatorial bases
    sage: sorted(A.basis_names())
    ['forest', 'fundamental_slide', 'glide', 'grothendieck', 'grove', 'key', 'lascoux', 'monomial', 'monomial_slide', 'schubert']
    sage: A.basis('lascoux')
    Polynomial ring in x0, x1, ... over Integer Ring in the Lascoux polynomial basis
    sage: A.elementary()
    Traceback (most recent call last):
    ...
    ValueError: the elementary symmetric basis depends on the number of variables: use PolynomialAlgebra(R, n)
    sage: TestSuite(A).run()

    sage: B = PolynomialAlgebra(ZZ, 2); B
    Polynomial ring in x0, x1 over Integer Ring with combinatorial bases
    sage: B.number_of_variables()
    2
    sage: 'elementary' in B.basis_names()
    True
    sage: TestSuite(B).run()

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.number_of_variables"></a>

#### number\_of\_variables

```python
def number_of_variables()
```

The number of variables, or ``None`` for the ring in infinitely many variables.

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.basis_names"></a>

#### basis\_names

```python
def basis_names()
```

The names of the bases available in this ring.

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.basis"></a>

#### basis

```python
@cached_method
def basis(name)
```

The realization called ``name`` (see :meth:`basis_names`).

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.monomial"></a>

#### monomial

```python
def monomial()
```

The monomial basis `x^\alpha`.

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: x = PolynomialAlgebra(QQ).monomial()
    sage: x[2, 0, 1] * x[0, 1]
    x[2, 1, 1]
    sage: x[2, 0, 1].expand()
    x0^2*x2

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.schubert"></a>

#### schubert

```python
def schubert()
```

Schubert polynomials `\mathfrak S_w`, indexed by permutations.

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: S = PolynomialAlgebra(QQ).schubert()
    sage: S(Permutation([3, 1, 2])) == S[3, 1, 2] == S([3, 1, 2, 4])
    True
    sage: S[3, 1, 2] * S[1, 3, 2]
    S[3, 2, 1] + S[4, 1, 2, 3]
    sage: S[2, 3, 1].expand()
    x0*x1
    sage: S[1, 3, 2].degree()
    1

In `n` variables the last descent must be at most `n`::

    sage: PolynomialAlgebra(QQ, 2).schubert()([1, 2, 4, 3])
    Traceback (most recent call last):
    ...
    ValueError: [1, 2, 4, 3] has a descent beyond the 2 variables of Polynomial ring in x0, x1 over Rational Field with combinatorial bases

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.grothendieck"></a>

#### grothendieck

```python
def grothendieck()
```

Grothendieck polynomials `\mathfrak G_w` (at `\beta = -1`, the classical convention), indexed by permutations.

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: A = PolynomialAlgebra(QQ); G = A.grothendieck()
    sage: G[1, 3, 2].expand()
    -x0*x1 + x0 + x1
    sage: A.schubert()(G[1, 3, 2])
    S[1, 3, 2] - S[2, 3, 1]
    sage: G[2, 1] * G[2, 1]
    G[3, 1, 2]
    sage: G[1, 3, 2] * G[2, 1]
    G[2, 3, 1] + G[3, 1, 2] - G[3, 2, 1]

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.key"></a>

#### key

```python
def key()
```

Key polynomials (Demazure characters) `\kappa_\alpha`.

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: k = PolynomialAlgebra(QQ).key()
    sage: k[0, 2].expand()
    x0^2 + x0*x1 + x1^2
    sage: k(KeyPolynomials(QQ)([0, 2]))
    k[0, 2]

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.fundamental_slide"></a>

#### fundamental\_slide

```python
def fundamental_slide()
```

Fundamental slide polynomials `\mathfrak F_\alpha` (Assaf-Searles).

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: A = PolynomialAlgebra(QQ); F = A.fundamental_slide()
    sage: F(A.schubert()[2, 1, 5, 3, 4])
    F[1, 0, 2] + F[2, 0, 1] + F[3]

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.monomial_slide"></a>

#### monomial\_slide

```python
def monomial_slide()
```

Monomial slide polynomials `\mathfrak M_\alpha` (Assaf-Searles).

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: A = PolynomialAlgebra(QQ); M = A.monomial_slide()
    sage: M(A.fundamental_slide()[1, 0, 2])
    M[1, 0, 2] + M[1, 1, 1]

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.forest"></a>

#### forest

```python
def forest()
```

Forest polynomials `\mathfrak P_F` (Nadeau-Spink-Tewari), indexed by weak compositions.

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: A = PolynomialAlgebra(QQ); P = A.forest()
    sage: P(A.schubert()[2, 1, 5, 3, 4])
    P[1, 0, 2]

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.glide"></a>

#### glide

```python
def glide()
```

Glide polynomials (Pechenik-Searles), the K-theoretic fundamental slides, at `\beta = -1`.

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: A = PolynomialAlgebra(QQ); Gl = A.glide()
    sage: Gl[0, 2].expand()
    -x0^2*x1 - x0*x1^2 + x0^2 + x0*x1 + x1^2
    sage: Gl(A.grothendieck()[1, 3, 2])
    Gl[0, 1]

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.lascoux"></a>

#### lascoux

```python
def lascoux()
```

Lascoux polynomials, the K-theoretic key polynomials, at `\beta = -1`.

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: A = PolynomialAlgebra(QQ); L = A.lascoux()
    sage: A.glide()(L[0, 2])
    Gl[0, 2]
    sage: L[1, 0, 2].expand()
    x0^2*x1^2*x2 + x0^2*x1*x2^2 - x0^2*x1^2 - 2*x0^2*x1*x2 - x0*x1^2*x2 - x0^2*x2^2 - x0*x1*x2^2 + x0^2*x1 + x0*x1^2 + x0^2*x2 + x0*x1*x2 + x0*x2^2

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.grove"></a>

#### grove

```python
def grove()
```

Grove polynomials, the K-theoretic forest polynomials, at `\beta = -1`.

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: A = PolynomialAlgebra(QQ); Gr = A.grove()
    sage: Gr[0, 2].expand()
    -x0^2*x1 - x0*x1^2 + x0^2 + x0*x1 + x1^2

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebra.elementary"></a>

#### elementary

```python
def elementary()
```

Products of elementary symmetric polynomials, a basis of the ring in `n` variables.

The index is a tuple `(a_1, \ldots, a_{n-1}, b_1, \ldots, b_r)`, `r \geq 1`, standing for
`\prod_{j=1}^{n-1} e_{a_j}(x_0, \ldots, x_{j-1}) \cdot \prod_i e_{b_i}(x_0, \ldots, x_{n-1})` with
`b_1 \leq \cdots \leq b_r` (a single `b_1 = 0` when there is no full-alphabet factor): the
basis `\{\prod_{j<n} e_{a_j}(x_0..x_{j-1})\}` of the polynomial ring as a free module over the
symmetric polynomials, times the monomials in `e_1, \ldots, e_n`. Since the meaning of an index
depends on `n`, this basis exists only in ``PolynomialAlgebra(R, n)``.

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: B = PolynomialAlgebra(QQ, 3); E = B.elementary()
    sage: E[1, 0, 2].expand()                 # e_1(x0) e_0(x0, x1) e_2(x0, x1, x2)
    x0^2*x1 + x0^2*x2 + x0*x1*x2
    sage: E[0, 0, 1, 1].expand()              # e_1(x0, x1, x2)^2
    x0^2 + 2*x0*x1 + x1^2 + 2*x0*x2 + 2*x1*x2 + x2^2
    sage: E([1, 0, 2, 0, 1]) == E[1, 0, 1, 2]  # e_0 factors drop out, full-alphabet degrees are sorted
    True
    sage: E(B.monomial()[2])
    -E[0, 2, 0] + E[1, 1, 0]
    sage: B.schubert()(E[0, 0, 1])
    S[1, 2, 4, 3]

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebraBasis"></a>

## PolynomialAlgebraBasis Objects

```python
class PolynomialAlgebraBasis(CombinatorialFreeModule)
```

A basis of :func:`PolynomialAlgebra`; subclasses fix the index set.

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: k = PolynomialAlgebra(QQ).key(); k
    Polynomial ring in x0, x1, ... over Rational Field in the key polynomial basis
    sage: k([2, 0, 1, 0, 0])
    k[2, 0, 1]
    sage: k.one()
    k[]
    sage: TestSuite(k).run()
    sage: TestSuite(PolynomialAlgebra(QQ, 2).schubert()).run()

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebraBasis.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(key)
```

``k[2, 0, 1]``, ``S[3, 1, 2]``: the basis element with that index.

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebraBasis.Element"></a>

## Element Objects

```python
class Element(CombinatorialFreeModule.Element)
```

<a id="schubmult.sage.polynomial_algebra.PolynomialAlgebraBasis.Element.expand"></a>

#### expand

```python
def expand()
```

The polynomial in ``x0, x1, ...``: in the `n` variables of the ring, or in as many as the
support needs (at least one).

EXAMPLES::

    sage: from schubmult.sage import PolynomialAlgebra
    sage: A = PolynomialAlgebra(QQ)
    sage: A.key()[1, 0, 2].expand()
    x0^2*x1 + x0*x1^2 + x0^2*x2 + x0*x1*x2 + x0*x2^2
    sage: A.lascoux()[0, 2].expand()
    -x0^2*x1 - x0*x1^2 + x0^2 + x0*x1 + x1^2
    sage: A.monomial().one().expand().parent()
    Multivariate Polynomial Ring in x0 over Rational Field
    sage: PolynomialAlgebra(QQ, 3).monomial().one().expand().parent()
    Multivariate Polynomial Ring in x0, x1, x2 over Rational Field

