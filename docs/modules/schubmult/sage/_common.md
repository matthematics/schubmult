<a id="schubmult.sage._common"></a>

# schubmult.sage.\_common

Shared machinery for Sage parents backed by schubmult rings.

A backed ring is a :class:`~sage.combinat.free_module.CombinatorialFreeModule` indexed by
permutations whose base ring is an infinite polynomial ring over the scalars in the coefficient
alphabets (``y``, ``z``, ``q``, ...). Arithmetic is delegated to a schubmult ring object
(:meth:`SchubmultBackedRing._schub_ring`) and coefficients are converted at the boundary.

Index conventions: Sage variables are 0-indexed (``x0``, ``y_0``, ``q_0``), schubmult's are 1-indexed
(``x_1``, ``y_1``, ``q_1``); see :mod:`schubmult.sage._convert`.

<a id="schubmult.sage._common.coefficient_into"></a>

#### coefficient\_into

```python
def coefficient_into(c, T, aliases=None)
```

Move a base-ring coefficient into the finite ring ``T`` (variables ``a<i>`` for ``a_<i>``, plus
named scalars like ``beta``); fractions land in ``T.fraction_field()``. ``aliases`` renames
base-ring variables (``{'beta_0': 'beta'}``).

<a id="schubmult.sage._common.SchubmultBackedElement"></a>

## SchubmultBackedElement Objects

```python
class SchubmultBackedElement(CombinatorialFreeModule.Element)
```

<a id="schubmult.sage._common.SchubmultBackedElement.project"></a>

#### project

```python
def project(n)
```

The image in the cohomology of the flag variety of `\CC^n` (equivariant, quantum, or partial as
the ring dictates): drop the terms indexed by permutations outside `S_n`, whose classes vanish there.

The rings are stable -- `\mathfrak S_w` is `\mathfrak S_w` for every `n` -- so a product
contains every class of the infinite flag variety; restricting to one `n` is a projection.

EXAMPLES::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: X = DoubleSchubertPolynomialRing(QQ)
    sage: f = X([2, 3, 1]) * X([3, 1, 2]); f
    (y_2-y_0)*X_y[3, 2, 1] + X_y[4, 2, 1, 3]
    sage: f.project(3)                # in H_T^*(Fl(3)) only S_{321} survives
    (y_2-y_0)*X_y[3, 2, 1]

<a id="schubmult.sage._common.SchubmultBackedElement.expand"></a>

#### expand

```python
def expand()
```

Expand into a polynomial in ``x0, x1, ...`` and the coefficient variables ``y0, q0, ...``.

There are `n` variables ``x``, `n` the size of the largest permutation involved (as for
:meth:`sage.combinat.schubert_polynomial.SchubertPolynomial_class.expand`); coefficient
letters get exactly the indices that occur.

<a id="schubmult.sage._common.SchubmultBackedElement.to_symmetric_function"></a>

#### to\_symmetric\_function

```python
def to_symmetric_function(n=None)
```

The symmetric function whose expansion in `x_0, \ldots, x_{n-1}` is this element.

``n`` defaults to the number of ``x`` variables that actually occur. Returns an element of
``SymmetricFunctions(C)`` in the Schur basis supported on partitions with at most ``n`` parts
(the unique such preimage), ``C`` the ring of the coefficient variables (``y``, ``q``, ``beta``,
...) that occur; raises ``ValueError`` if the expansion is not symmetric in those variables.
Grassmannian Schubert polynomials with descent at `n` are the Schur functions
`s_\lambda(x_0, \ldots, x_{n-1})`; double ones are factorial Schur functions.

EXAMPLES::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing, GrothendieckPolynomialRing
    sage: X = DoubleSchubertPolynomialRing(QQ)
    sage: X([2, 4, 1, 3]).to_symmetric_function()
    -(y0^2*y1+y0^2*y2)*s[] + (y0^2+y0*y1+y0*y2)*s[1] - (y0+y1+y2)*s[1, 1] - y0*s[2] + s[2, 1]
    sage: X([3, 2, 1]).to_symmetric_function()
    Traceback (most recent call last):
    ...
    ValueError: X_y[3, 2, 1] is not symmetric in x0, x1
    sage: X([3, 1, 2]).to_symmetric_function()  # a polynomial in x0 alone
    y0*y1*s[] - (y0+y1)*s[1] + s[2]
    sage: X([3, 1, 2]).to_symmetric_function(2)
    Traceback (most recent call last):
    ...
    ValueError: X_y[3, 1, 2] is not symmetric in x0, x1
    sage: G = GrothendieckPolynomialRing(ZZ)
    sage: G([1, 3, 2]).to_symmetric_function()
    s[1] + beta*s[1, 1]

<a id="schubmult.sage._common.SchubmultBackedRing"></a>

## SchubmultBackedRing Objects

```python
class SchubmultBackedRing(CombinatorialFreeModule)
```

Base class: subclasses set ``_alphabet`` (basis alphabet letter or ``None``), ``_alphabets`` (all
coefficient letters, sorted), and implement ``_schub_ring(alphabet)`` and ``_check_basis_perm``.

<a id="schubmult.sage._common.SchubmultBackedRing.from_symmetric_function"></a>

#### from\_symmetric\_function

```python
def from_symmetric_function(f, n)
```

Expand the symmetric function ``f`` in ``n`` variables `x_0, \ldots, x_{n-1}` in this basis.

A Schur function `s_\lambda(x_0, \ldots, x_{n-1})` is the Schubert polynomial of the Grassmannian
permutation with descent at `n` and shape `\lambda`; in the double ring the same polynomial
expands with coefficients in the second alphabet.

EXAMPLES::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: X = DoubleSchubertPolynomialRing(QQ); s = SymmetricFunctions(QQ).s()
    sage: X.from_symmetric_function(s[2, 1], 2)
    (y_1^2*y_0+y_1*y_0^2)*X_y[1] + (y_2*y_0+y_1*y_0+y_0^2)*X_y[1, 3, 2] + y_0*X_y[1, 4, 2, 3] + (y_2+y_1+y_0)*X_y[2, 3, 1] + X_y[2, 4, 1, 3]
    sage: SchubertPolynomialRing(QQ)(s[2, 1].expand(2, alphabet=['x0', 'x1']))
    X[2, 4, 1, 3]
    sage: X.from_symmetric_function(s[2, 1], 2).to_symmetric_function()
    s[2, 1]

