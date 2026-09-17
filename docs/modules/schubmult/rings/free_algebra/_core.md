<a id="schubmult.rings.free_algebra._core"></a>

# schubmult.rings.free\_algebra.\_core

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement"></a>

## FreeAlgebraElement Objects

```python
class FreeAlgebraElement(BaseRingElement)
```

Element of a free algebra, stored as a dict mapping basis keys to coefficients.

Keys are tuples of nonnegative integers (words in the word basis) or
basis-specific keys depending on the parent ring's basis. Supports
arithmetic operations, basis changes, and word-level operations like
injection, prefix, suffix, and interval extraction.

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

Compute the inner product of this element with a polynomial.

Converts to WordBasis and pairs coefficient-by-coefficient with
the monomial expansion of *poly* in *genset*.

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

Compute the pairing of this element with *other* via the monomial basis.

Converts *self* to WordBasis and *other* to MonomialBasis, then
sums products of matching coefficients.

**Arguments**:

- `other` - Another element to pair with.
  

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

Expand a word tuple into the double Schubert basis via Pieri products.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.tup_expand"></a>

#### tup\_expand

```python
@staticmethod
@cache
def tup_expand(tup)
```

Expand a word tuple into the single Schubert basis via divide-and-conquer Pieri products.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.change_basis"></a>

#### change\_basis

```python
def change_basis(other_basis)
```

Convert this element to another basis.

**Arguments**:

- `other_basis` - The target basis class (e.g. WordBasis, SchubertBasis).
  

**Returns**:

  A new FreeAlgebraElement in the target basis's ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.schub_expand"></a>

#### schub\_expand

```python
def schub_expand()
```

Expand this element into a single Schubert polynomial ring element.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.schub_double_expand"></a>

#### schub\_double\_expand

```python
def schub_double_expand()
```

Expand this element into a double Schubert polynomial ring element.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.bcoproduct"></a>

#### bcoproduct

```python
def bcoproduct()
```

Compute the bar-coproduct of this element in the tensor ring.

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

Free algebra ring with a configurable basis.

The algebra operates on :class:`FreeAlgebraElement` instances whose keys
are determined by the chosen basis (default :class:`WordBasis`). Supports
multiplication, tensor products, coproducts, and basis changes.

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

Internal product (``@`` operator) or scalar multiplication.

If *other* is a scalar, multiplies all coefficients. If *other* is
a FreeAlgebraElement, computes the internal product via the basis.

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

Skew schubert by elem sym

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

