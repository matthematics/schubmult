<a id="schubmult.sage.quantum_schubert"></a>

# schubmult.sage.quantum\_schubert

Quantum and quantum double Schubert polynomials

The quantum Schubert polynomials `\mathfrak{S}^q_w(x)` of Fomin-Gelfand-Postnikov represent Schubert
classes in the (small) quantum cohomology ring `QH^*(Fl_n)`; the quantum double Schubert polynomials
`\mathfrak{S}^q_w(x; y)` (Kirillov-Maeno, Ciocan-Fontanine-Fulton) do the same equivariantly. With
block sizes ``parabolic = (n_1, ..., n_k)`` one gets the quantum cohomology of the partial flag
variety `Fl(n_1, n_1 + n_2, \ldots)`: the basis is indexed by permutations whose descents lie at the
block boundaries `N_j = n_1 + \cdots + n_j`, and products are computed in the full flag ring and
projected by the Peterson-Woodward comparison formula.

The recorded blocks are followed by an implicit last block that is never recorded and is as large as
any computation needs it to be. Concretely, a basis permutation may have descents exactly at the
recorded boundaries `N_1, \ldots, N_k` (so `N_k` itself is allowed) and must be increasing beyond
`N_k`; a product whose result reaches past `N_k` is computed with one extra block appended, and the
answer does not depend on how large that block is made. The recorded block sizes are never changed:
enlarging the last recorded block would change which `q` has which degree, and adding singleton blocks
instead gives a different (and, as polynomials, inconsistent) ring.

For blocks `(2, 3)` (the Grassmannian `Gr(2, 5)`) the basis elements inside the `2 \times 3` box are
literally Schur polynomials in `x_0, x_1`, and `q` enters a product only through classes that need the
implicit block; those carry `q`'s in their polynomial that cancel again on expansion::

    sage: from schubmult.sage import QuantumSchubertPolynomialRing
    sage: G = QuantumSchubertPolynomialRing(QQ, parabolic=(2, 3))
    sage: s32, s1 = G([3, 5, 1, 2, 4]), G([1, 3, 2])      # sigma_32, sigma_1
    sage: s32.expand(), s1.expand()
    (x0^3*x1^2 + x0^2*x1^3, x0 + x1)
    sage: s32 * s1                                        # sigma_33 + q sigma_1  (+ sigma_42, which lives past the box)
    q_0*Xq[1, 3, 2] + Xq[3, 6, 1, 2, 4, 5] + Xq[4, 5, 1, 2, 3]
    sage: G([3, 6, 1, 2, 4, 5]).expand()                  # sigma_42 = s_42 - q_0 s_1
    x0^4*x1^2 + x0^3*x1^3 + x0^2*x1^4 - x0*q0 - x1*q0
    sage: (s32 * s1).expand() == s32.expand() * s1.expand()
    True

The quantum parameters are ``q_0, q_1, ...`` (0-indexed like everything on the Sage side; the
schubmult parameter `q_i` is ``q_{i-1}``); the base ring is ``R[q]`` resp. ``R[q, y, z]``.

EXAMPLES::

    sage: from schubmult.sage import QuantumSchubertPolynomialRing, QuantumDoubleSchubertPolynomialRing
    sage: Q = QuantumSchubertPolynomialRing(QQ); Q
    Quantum Schubert polynomial ring with Xq basis over Rational Field
    sage: Q([2, 1]) * Q([2, 1])
    q_0*Xq[1] + Xq[3, 1, 2]
    sage: Q([3, 1, 2]).expand()
    x0^2 - q0

The quantum Monk formula `\mathfrak{S}^q_{s_1} \mathfrak{S}^q_{s_1} = \mathfrak{S}^q_{312} + q_1`
is visible above; at `q = 0` the classical product returns::

    sage: (Q([2, 1]) * Q([2, 1])).expand().subs(q0=0)
    x0^2

Products agree with polynomial multiplication in `\QQ[q][x]`::

    sage: f = Q([3, 1, 2]) * Q([2, 3, 1])
    sage: f.expand() == Q([3, 1, 2]).expand() * Q([2, 3, 1]).expand()
    True

Quantum double::

    sage: QD = QuantumDoubleSchubertPolynomialRing(QQ)
    sage: QD([2, 1]) * QD([2, 1])
    q_0*Xq_y[1] + (y_1-y_0)*Xq_y[2, 1] + Xq_y[3, 1, 2]
    sage: QD([3, 1, 2]).expand()
    x0^2 - x0*y0 - x0*y1 + y0*y1 - q0

Parabolic (block sizes 2, 3, so the flag variety `Fl(2, 5)`; the position-5 boundary and anything
increasing beyond it are allowed)::

    sage: P = QuantumSchubertPolynomialRing(QQ, parabolic=(2, 3)); P
    Parabolic quantum Schubert polynomial ring with Xq basis for block sizes (2, 3) over Rational Field
    sage: P([2, 1, 3])
    Traceback (most recent call last):
    ...
    ValueError: [2, 1] is not parabolic for block sizes (2, 3): descents must lie in [2, 5]
    sage: P([2, 3, 1]) * P([2, 4, 1, 3])
    Xq[3, 5, 1, 2, 4]
    sage: P([1, 3, 2, 4, 6, 5])  # descent at 5 = N_2, increasing afterwards
    Xq[1, 3, 2, 4, 6, 5]

<a id="schubmult.sage.quantum_schubert.QuantumSchubertPolynomialRing"></a>

#### QuantumSchubertPolynomialRing

```python
def QuantumSchubertPolynomialRing(R, parabolic=None)
```

Return the ring of quantum Schubert polynomials `\mathfrak{S}^q_w(x)` over ``R``.

INPUT:

- ``R`` -- a commutative ring; the base ring of the result is ``R[q_0, q_1, ...]``
- ``parabolic`` -- (optional) block sizes `(n_1, \ldots, n_k)` of a partial flag variety; basis
  permutations may only have descents at the block boundaries `n_1 + \cdots + n_j` (an implicit
  unbounded last block follows the recorded ones)

EXAMPLES::

    sage: from schubmult.sage import QuantumSchubertPolynomialRing
    sage: Q = QuantumSchubertPolynomialRing(ZZ); Q
    Quantum Schubert polynomial ring with Xq basis over Integer Ring
    sage: Q.base_ring()
    Infinite polynomial ring in q over Integer Ring
    sage: TestSuite(Q).run()
    sage: Q([1, 3, 2]) * Q([2, 1])
    Xq[2, 3, 1] + Xq[3, 1, 2]
    sage: f = Q([2, 3, 1]) * Q([3, 1, 2]); f
    q_1*q_0*Xq[1] + Xq[4, 2, 1, 3]
    sage: f.expand() == Q([2, 3, 1]).expand() * Q([3, 1, 2]).expand()
    True

Ordinary Schubert polynomials coerce in (as polynomials in `x`, re-expanded in the quantum basis)::

    sage: Q(SchubertPolynomialRing(ZZ)([3, 1, 2]))
    q_0*Xq[1] + Xq[3, 1, 2]

Parabolic, blocks `(1, 2)`: the projective plane `\PP^2 = Gr(1, 3)` with the implicit block after it.
`\sigma_1^2 = \sigma_2` has no `q`; `\sigma_1 \sigma_2 = q` in `QH^*(\PP^2)`, and the class
`\sigma_3` that only exists past the box appears alongside (its polynomial is `x_0^3 - q_0`)::

    sage: P = QuantumSchubertPolynomialRing(QQ, parabolic=(1, 2))
    sage: TestSuite(P).run()
    sage: P([2, 1, 3]) * P([2, 1, 3])
    Xq[3, 1, 2]
    sage: P([2, 1, 3]) * P([3, 1, 2])
    q_0*Xq[1] + Xq[4, 1, 2, 3]
    sage: P([4, 1, 2, 3]).expand()
    x0^3 - q0

<a id="schubmult.sage.quantum_schubert.QuantumDoubleSchubertPolynomialRing"></a>

#### QuantumDoubleSchubertPolynomialRing

```python
def QuantumDoubleSchubertPolynomialRing(R,
                                        alphabet="y",
                                        coefficient_alphabets=("y", "z"),
                                        parabolic=None)
```

Return the ring of quantum double Schubert polynomials `\mathfrak{S}^q_w(x; \text{alphabet})` over ``R``.

INPUT:

- ``R`` -- a commutative ring; the base ring of the result is ``R[q, alphabets]``
- ``alphabet`` -- (default: ``'y'``) the letter of the second alphabet of the basis elements
- ``coefficient_alphabets`` -- (default: ``('y', 'z')``) letters available in coefficients
- ``parabolic`` -- (optional) block sizes of a partial flag variety, as for
  :func:`QuantumSchubertPolynomialRing`

EXAMPLES::

    sage: from schubmult.sage import QuantumDoubleSchubertPolynomialRing
    sage: QD = QuantumDoubleSchubertPolynomialRing(QQ); QD
    Quantum double Schubert polynomial ring in the alphabet y with Xq_y basis over Rational Field
    sage: QD.base_ring()
    Infinite polynomial ring in q, y, z over Rational Field
    sage: TestSuite(QD).run()
    sage: QD([3, 1, 2]) * QD([2, 1])
    q_0*Xq_y[1, 3, 2] + (y_2-y_0)*Xq_y[3, 1, 2] + Xq_y[4, 1, 2, 3]

Mixed alphabets, as for :func:`~schubmult.sage.DoubleSchubertPolynomialRing`::

    sage: QZ = QuantumDoubleSchubertPolynomialRing(QQ, 'z')
    sage: QD([2, 1]) * QZ([2, 1])
    q_0*Xq_y[1] + (y_1-z_0)*Xq_y[2, 1] + Xq_y[3, 1, 2]

Products agree with polynomial multiplication::

    sage: f = QD([3, 1, 2]) * QD([2, 3, 1])
    sage: f.expand() == QD([3, 1, 2]).expand() * QD([2, 3, 1]).expand()
    True

Parabolic quantum double, blocks `(1, 2)`::

    sage: PD = QuantumDoubleSchubertPolynomialRing(QQ, parabolic=(1, 2))
    sage: TestSuite(PD).run()
    sage: f = PD([2, 1, 3]) * PD([2, 1, 3]); f
    (y_1-y_0)*Xq_y[2, 1] + Xq_y[3, 1, 2]
    sage: f.expand() == PD([2, 1, 3]).expand()^2
    True

<a id="schubmult.sage.quantum_schubert.QuantumSchubertPolynomial_class"></a>

## QuantumSchubertPolynomial\_class Objects

```python
class QuantumSchubertPolynomial_class(SchubmultBackedElement)
```

<a id="schubmult.sage.quantum_schubert.QuantumSchubertPolynomial_class.expand"></a>

#### expand

```python
def expand()
```

Expand into a polynomial in ``x0, x1, ...`` and ``q0, q1, ...``.

EXAMPLES::

    sage: from schubmult.sage import QuantumSchubertPolynomialRing
    sage: Q = QuantumSchubertPolynomialRing(ZZ)
    sage: [Q(p).expand() for p in Permutations(3)]
    [1, x0 + x1, x0, x0*x1 + q0, x0^2 - q0, x0^2*x1 + x0*q0]
    sage: Q([2, 1]).expand().parent()
    Multivariate Polynomial Ring in x0, x1 over Integer Ring

<a id="schubmult.sage.quantum_schubert._QuantumMixin"></a>

## \_QuantumMixin Objects

```python
class _QuantumMixin()
```

<a id="schubmult.sage.quantum_schubert._QuantumMixin.parabolic"></a>

#### parabolic

```python
def parabolic()
```

The block sizes of the partial flag variety, or ``None`` for the full flag variety.

EXAMPLES::

    sage: from schubmult.sage import QuantumSchubertPolynomialRing
    sage: QuantumSchubertPolynomialRing(QQ).parabolic() is None
    True
    sage: QuantumSchubertPolynomialRing(QQ, parabolic=[2, 3]).parabolic()
    (2, 3)

<a id="schubmult.sage.quantum_schubert.QuantumSchubertPolynomialRing_xbasis"></a>

## QuantumSchubertPolynomialRing\_xbasis Objects

```python
class QuantumSchubertPolynomialRing_xbasis(_QuantumMixin, SchubmultBackedRing)
```

<a id="schubmult.sage.quantum_schubert.QuantumSchubertPolynomialRing_xbasis.__init__"></a>

#### \_\_init\_\_

```python
def __init__(R, parabolic)
```

EXAMPLES::

    sage: from schubmult.sage import QuantumSchubertPolynomialRing
    sage: Q = QuantumSchubertPolynomialRing(QQ, parabolic=(2, 3))
    sage: Q == loads(dumps(Q))
    True
    sage: Q is QuantumSchubertPolynomialRing(QQ, parabolic=[2, 3])
    True

<a id="schubmult.sage.quantum_schubert.QuantumSchubertPolynomialRing_xbasis.product_on_basis"></a>

#### product\_on\_basis

```python
def product_on_basis(left, right)
```

EXAMPLES::

    sage: from schubmult.sage import QuantumSchubertPolynomialRing
    sage: Q = QuantumSchubertPolynomialRing(QQ)
    sage: Q.product_on_basis(Permutation([3, 1, 2]), Permutation([2, 1]))
    q_0*Xq[1, 3, 2] + Xq[4, 1, 2, 3]

<a id="schubmult.sage.quantum_schubert.QuantumDoubleSchubertPolynomialRing_xbasis"></a>

## QuantumDoubleSchubertPolynomialRing\_xbasis Objects

```python
class QuantumDoubleSchubertPolynomialRing_xbasis(_QuantumMixin,
                                                 SchubmultBackedRing)
```

<a id="schubmult.sage.quantum_schubert.QuantumDoubleSchubertPolynomialRing_xbasis.__init__"></a>

#### \_\_init\_\_

```python
def __init__(R, alphabet, alphabets, parabolic)
```

EXAMPLES::

    sage: from schubmult.sage import QuantumDoubleSchubertPolynomialRing
    sage: QD = QuantumDoubleSchubertPolynomialRing(QQ, 'z', parabolic=(1, 2))
    sage: QD == loads(dumps(QD))
    True

<a id="schubmult.sage.quantum_schubert.QuantumDoubleSchubertPolynomialRing_xbasis.alphabet"></a>

#### alphabet

```python
def alphabet()
```

The letter of the second alphabet of the basis elements.

EXAMPLES::

    sage: from schubmult.sage import QuantumDoubleSchubertPolynomialRing
    sage: QuantumDoubleSchubertPolynomialRing(QQ, 'z').alphabet()
    'z'

