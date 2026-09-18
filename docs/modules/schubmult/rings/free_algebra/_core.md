<a id="schubmult.rings.free_algebra._core"></a>

# schubmult.rings.free\_algebra.\_core

`FreeAlgebra` and `FreeAlgebraElement`: the graded dual of the polynomial algebra.

See the package docstring (`schubmult.rings.free_algebra`) for the duality. This module
holds the ring and element classes; the individual bases live in sibling modules and
plug in through the `FreeAlgebraBasis` interface. ``FA``, ``ASx``, ``AGx`` are the
standard instances.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement"></a>

## FreeAlgebraElement Objects

```python
class FreeAlgebraElement(BaseRingElement)
```

An element of a `FreeAlgebra`: ``{basis_key: coefficient}``.

In the `WordBasis` a key is a word -- a tuple of nonnegative integers -- dual to
the monomial with that exponent vector. In other bases the key is that basis's
combinatorial index together with a number of variables (e.g. ``(perm, numvars)``
for `SchubertBasis`). Beyond ring arithmetic, elements support basis changes
(``change_basis``), the duality pairing with polynomials (``pairing``,
``poly_inner_product``), expansion into Schubert rings (``schub_expand``), and
word-level operations (``inject``, ``prefix``, ``suffix``, ``interval``, ``split``,
``factorize``) that are computed in the word basis and transported back.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.interleave"></a>

#### interleave

```python
def interleave(other, zero_pad=True)
```

Interleave two elements letter-by-letter in the word basis.

Converts both elements to WordBasis, then interleaves each pair of
words by alternating entries (a1, b1, a2, b2, ...). Shorter words
are zero-padded when ``zero_pad`` is True.

**Arguments**:

- `other` - Another FreeAlgebraElement to interleave with.
- `zero_pad` - If True, pad shorter words with zeros.
  

**Returns**:

  The interleaved element in the original basis.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.inject"></a>

#### inject

```python
def inject(i, other)
```

Insert another element's words at position *i* in this element's words.

Delegates to the current basis's ``inject`` classmethod.

**Arguments**:

- `i` - Nonnegative integer insertion index.
- `other` - Another FreeAlgebraElement to inject.
  

**Returns**:

  A new element with the injected words.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.prefix"></a>

#### prefix

```python
def prefix(length)
```

Extract the first *length* letters of each word.

Delegates to the current basis's ``prefix`` classmethod.

**Arguments**:

- `length` - Nonnegative integer prefix length.
  

**Returns**:

  A new element containing the prefixes.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.suffix"></a>

#### suffix

```python
def suffix(length)
```

Extract the last *length* letters of each word.

Delegates to the current basis's ``suffix`` classmethod.

**Arguments**:

- `length` - Nonnegative integer suffix length.
  

**Returns**:

  A new element containing the suffixes.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.interval"></a>

#### interval

```python
def interval(start, stop)
```

Extract a subword from position *start* to *stop* in each word.

Delegates to the current basis's ``interval`` classmethod.

**Arguments**:

- `start` - Nonnegative start index (inclusive).
- `stop` - Stop index (exclusive).
  

**Returns**:

  A new element containing the subwords.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.poly_inner_product"></a>

#### poly\_inner\_product

```python
def poly_inner_product(poly, genset, n)
```

The duality pairing of this element with a polynomial expression.

Expands ``poly`` into monomials in ``genset`` (exponent vectors padded/truncated to
``n`` variables, or trailing zeros stripped if ``n`` is ``None``), converts ``self``
to the word basis, and sums ``coeff_word * coeff_monomial`` over matching
word/exponent-vector pairs.

**Arguments**:

- `poly` - A polynomial expression.
- `genset` - The generating set of variables for the polynomial.
- `n` - Number of variables to use (or None for automatic).
  

**Returns**:

  The integer inner product value.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.kill_zero"></a>

#### kill\_zero

```python
def kill_zero(fat=False, val=S.Zero)
```

Remove zeros from each word key.

In the word basis, strips all zero entries from each key. When
*fat* is True, multiplies the coefficient by ``val`` raised to
the number of removed zeros instead of simply dropping them.

**Arguments**:

- `fat` - If True, weight by ``val`` per removed zero.
- `val` - The value to raise per zero when *fat* is True.
  

**Returns**:

  A new element with zeros removed, in the original basis.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.as_coefficients_dict"></a>

#### as\_coefficients\_dict

```python
def as_coefficients_dict()
```

Return a dict mapping printing terms to sympified coefficients.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.expand"></a>

#### expand

```python
def expand(deep=True, *args, **kwargs)
```

Expand all coefficients symbolically.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.pairing"></a>

#### pairing

```python
def pairing(other)
```

The duality pairing with a `PolynomialAlgebraElement`.

Converts ``self`` to the word basis and ``other`` to the monomial basis, then sums
``coeff_word * coeff_monomial`` over words equal to exponent vectors. This is the
pairing under which the free algebra is the graded dual of the polynomial algebra.

**Arguments**:

- `other` - A polynomial algebra element to pair with.
  

**Returns**:

  The integer pairing value.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.hom_nsym"></a>

#### hom\_nsym

```python
def hom_nsym()
```

Apply the homomorphism to noncommutative symmetric functions.

Strips zero entries from each word key (dropping them from the word)
and returns the result in the original basis.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.tup_double_expand"></a>

#### tup\_double\_expand

```python
@staticmethod
@cache
def tup_double_expand(tup)
```

Realize the word ``tup`` as a double Schubert (separated-descents) ring element: the
product ``prod_i S_{uncode([tup[i]])}`` with each factor in its own variable.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.tup_expand"></a>

#### tup\_expand

```python
@staticmethod
@cache
def tup_expand(tup)
```

Realize the word ``tup`` as a single Schubert (separated-descents) ring element: the
product ``prod_i S_{uncode([tup[i]])}`` with each factor in its own variable, computed
by divide-and-conquer. This is the map word -> ``h_{a_1}(x_1) h_{a_2}(x_2) ...`` that
sends the word basis onto complete-symmetric-in-one-variable products.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.change_basis"></a>

#### change\_basis

```python
def change_basis(other_basis)
```

Re-express this element in another `FreeAlgebraBasis`.

Uses ``self.ring._basis.transition(other_basis)``, which most bases implement by
routing through the `WordBasis`.

**Arguments**:

- `other_basis` - The target basis class (e.g. WordBasis, SchubertBasis).
  

**Returns**:

  A new FreeAlgebraElement in the target basis's ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.schub_expand"></a>

#### schub\_expand

```python
def schub_expand()
```

Realize this element in the single Schubert separated-descents ring via ``tup_expand``
(each word becomes a product of one-variable complete symmetric functions).

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.schub_double_expand"></a>

#### schub\_double\_expand

```python
def schub_double_expand()
```

Double-alphabet analogue of ``schub_expand`` via ``tup_double_expand``.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.bcoproduct"></a>

#### bcoproduct

```python
def bcoproduct()
```

The "bar" coproduct (see `WordBasis.bcoproduct`) in the tensor square ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.factorize"></a>

#### factorize

```python
def factorize(j)
```

Split each word at position *j*, returning a tensor element.

In the word basis, each word ``w`` maps to ``(w[:j], w[j:])``.
The result is expressed in the tensor ring of the original basis.

**Arguments**:

- `j` - Position at which to split each word.
  

**Returns**:

  An element of the tensor product ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.remove_zeros"></a>

#### remove\_zeros

```python
def remove_zeros(inserter=S.One)
```

Remove zero entries from each key, weighting by *inserter* per zero removed.

**Arguments**:

- `inserter` - Scalar multiplied per removed zero (default 1).
  

**Returns**:

  A new element with zeros stripped from keys.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.__truediv__"></a>

#### \_\_truediv\_\_

```python
def __truediv__(other)
```

Skew by a permutation: ``elem / u`` applies ``skew_element(w, u, n)`` to each ``(w, n)`` key
(the dual of multiplying by ``S_u`` on the polynomial side).

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.split"></a>

#### split

```python
def split(p)
```

Split each word at position *p* into a tensor element.

Words shorter than *p* are placed entirely in the left factor.

**Arguments**:

- `p` - Position at which to split.
  

**Returns**:

  An element of the tensor product ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.to_schub"></a>

#### to\_schub

```python
def to_schub(sym=False)
```

Convert to a Schubert ring element via ``tup_to_schub``.

**Arguments**:

- `sym` - If True, use symmetric expansion.
  

**Returns**:

  An element of the single Schubert ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra"></a>

## FreeAlgebra Objects

```python
class FreeAlgebra(BaseRing)
```

The free algebra on generators indexed by nonnegative integers, in a chosen basis.

The ring is basis-agnostic; a `FreeAlgebraBasis` *class* (not instance) supplies the
key type, product, coproduct, and transitions. ``FreeAlgebra(WordBasis)`` is the
concatenation algebra on words; ``FreeAlgebra(SchubertBasis)`` is the same algebra
written in the basis dual to Schubert polynomials. See the package docstring for
the duality with the polynomial algebra.

**Arguments**:

- `basis` - The basis class to use (default ``WordBasis``).
- `domain` - Coefficient domain (default ``EXRAW``).

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply every coefficient of *elem* by the scalar *x*.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.tensor_schub_expand"></a>

#### tensor\_schub\_expand

```python
def tensor_schub_expand(tensor)
```

Expand a tensor element into the Schubert polynomial tensor ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.__init__"></a>

#### \_\_init\_\_

```python
def __init__(basis=WordBasis, domain=None)
```

Initialize a FreeAlgebra with the given basis and coefficient domain.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.right_pad"></a>

#### right\_pad

```python
@staticmethod
def right_pad(tup, n)
```

Right-pad *tup* with zeros to length *n*.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
@cache
def coproduct_on_basis(key)
```

Compute the coproduct of a single basis key in the tensor ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.bcoproduct_on_basis"></a>

#### bcoproduct\_on\_basis

```python
@cache
def bcoproduct_on_basis(key)
```

Compute the bar-coproduct of a single basis key in the tensor ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.mul"></a>

#### mul

```python
def mul(elem, other)
```

Multiply two elements via the basis product rule.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.from_rc_graph"></a>

#### from\_rc\_graph

```python
def from_rc_graph(rc_graph)
```

Create an element from an RC graph.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.matmul"></a>

#### matmul

```python
def matmul(elem, other)
```

Internal (Kronecker) product (``@`` operator) or scalar multiplication.

If ``other`` is a scalar, multiplies all coefficients. If it is a `FreeAlgebraElement`,
computes the basis's ``internal_product`` -- essentially the Kronecker product of
noncommutative symmetric functions, enumerated in the word basis by nonnegative
integer matrices with prescribed row/column sums (see `WordBasis.internal_product`;
requires SageMath).

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.new"></a>

#### new

```python
def new(*x)
```

Create a new element from the given key or arguments.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Return the display symbol for basis key *k*.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.from_dict"></a>

#### from\_dict

```python
def from_dict(element)
```

Construct an element from a dict of ``{key: coefficient}`` pairs.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.skew_element"></a>

#### skew\_element

```python
def skew_element(w, u, n)
```

The skew element ``S_w / S_u`` in ``n`` variables (dual to multiplication by ``S_u``); see `SchubertBasis.skew_element`.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce a raw value into the coefficient domain.

<a id="schubmult.rings.free_algebra._core.make_AGx"></a>

#### make\_AGx

```python
def make_AGx(beta=None)
```

Create a Grothendieck-basis free algebra, optionally with custom beta.

