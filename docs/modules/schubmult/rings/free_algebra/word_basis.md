<a id="schubmult.rings.free_algebra.word_basis"></a>

# schubmult.rings.free\_algebra.word\_basis

`WordBasis`: the word (concatenation) basis of the free algebra, dual to the monomial basis.

A key is a word ``(a_1, ..., a_n)`` of nonnegative integers, dual to the monomial
``x_1^{a_1} ... x_n^{a_n}``. The product is concatenation; the coproduct splits each
letter ``a`` into ``(i, a - i)`` (dual to polynomial multiplication). This is the hub
basis: every other `FreeAlgebraBasis` implements its operations by transitioning to
words and back, and this module holds the ``transition_*`` routines from words into
each of the other bases.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis"></a>

## WordBasis Objects

```python
class WordBasis(FreeAlgebraBasis)
```

Word basis of the free algebra: keys are tuples of nonnegative integers (words), each
dual to the monomial whose exponent vector is that word. Product is concatenation.
See the module docstring.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Return the length vector of the RC graph as a word key.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Concatenate two words (dual to the variable-splitting coproduct on polynomials).

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.inject"></a>

#### inject

```python
@classmethod
def inject(cls, key1, i, key2, coeff=S.One)
```

Insert *key2* into *key1* at position *i*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.prefix"></a>

#### prefix

```python
@classmethod
def prefix(cls, key, length, coeff=S.One)
```

Return the first *length* letters of *key*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.suffix"></a>

#### suffix

```python
@classmethod
def suffix(cls, key, length, coeff=S.One)
```

Return the last *length* letters of *key*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.interval"></a>

#### interval

```python
@classmethod
def interval(cls, key, start, stop, coeff=S.One)
```

Return the subword ``key[start:stop]``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key, coeff=S.One)
```

The coproduct of a word, dual to polynomial multiplication.

Each letter ``a`` splits into all ``(i, a - i)`` pairs (``x_j^a`` is the sum over
ways to write it as ``x_j^i * x_j^{a-i}``); the pieces are combined letterwise
by a divide-and-conquer tensor product.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.bcoproduct"></a>

#### bcoproduct

```python
@classmethod
@cache
def bcoproduct(cls, key, coeff=S.One)
```

The "bar" coproduct: like `coproduct` but a zero letter is dropped rather than kept as
a ``0`` in the word, so word lengths are not preserved.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.try_internal_product"></a>

#### try\_internal\_product

```python
@classmethod
def try_internal_product(cls, key1, key2, coeff=S.One)
```

Compute the internal product via integer matrices (requires SageMath).

Uses shifted keys (incremented by 1) with ``IntegerMatrices``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.internal_product"></a>

#### internal\_product

```python
@classmethod
def internal_product(cls, key1, key2, coeff=S.One)
```

The internal (Kronecker) product of two compositions (words without zeros), as in
noncommutative symmetric functions: sum over nonnegative integer matrices with row
sums ``key1`` and column sums ``key2`` of the word read off the nonzero entries.
Requires SageMath's ``IntegerMatrices``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a bracket-notation symbol like ``[210]`` for the word *k*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.tup_expand"></a>

#### tup\_expand

```python
@staticmethod
@cache
def tup_expand(tup)
```

Expand a word tuple into the single Schubert basis via divide-and-conquer.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.jbasis_tup_expand"></a>

#### jbasis\_tup\_expand

```python
@staticmethod
@cache
def jbasis_tup_expand(tup)
```

Expand a word tuple into the Z basis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a word key to the Schubert basis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_jbasis"></a>

#### transition\_jbasis

```python
@classmethod
def transition_jbasis(cls, key)
```

Transition a word key to the J basis via Pieri products.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_jtbasis"></a>

#### transition\_jtbasis

```python
@classmethod
def transition_jtbasis(cls, key)
```

Transition a word key to the JT basis via normalization.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_forest"></a>

#### transition\_forest

```python
@classmethod
def transition_forest(cls, key)
```

Transition a word key to the forest basis via RC graph enumeration.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_grove"></a>

#### transition\_grove

```python
@classmethod
def transition_grove(cls, key)
```

Transition a word key to the grove basis via WC graph enumeration.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the MonomialBasis as the dual of WordBasis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_monomial_slide"></a>

#### transition\_monomial\_slide

```python
@classmethod
def transition_monomial_slide(cls, key)
```

Transition a word key to the monomial slide basis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_zbasis"></a>

#### transition\_zbasis

```python
@classmethod
def transition_zbasis(cls, key)
```

Expand a word in `ZBasis` by triangular elimination: repeatedly peel off the smallest
remaining word, subtracting the word expansion of the corresponding Z element.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_nelementary"></a>

#### transition\_nelementary

```python
@classmethod
def transition_nelementary(cls, tup)
```

Expand a composition in `NElementaryBasis`: signed sum over its refinements
(``(-1)^(|tup| - len(beta))``), via SageMath's ``Composition.finer``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_key"></a>

#### transition\_key

```python
@classmethod
def transition_key(cls, key)
```

Expand a word in `KeyBasis`: count RC graphs of length vector ``key`` whose extremal
weight equals their permutation's padded code (the dual of the key-to-monomial expansion).

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_lascoux"></a>

#### transition\_lascoux

```python
@classmethod
def transition_lascoux(cls, key)
```

Expand a word in `LascouxBasis`: the K-theoretic analogue of `transition_key` using WC graphs.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_glide"></a>

#### transition\_glide

```python
@classmethod
def transition_glide(cls, key)
```

Expand a word in `GlideBasis`: for each WC graph of weight ``key``, take the length vector of
its ``dst`` (destandardization); the first graph seen at each weight is the representative
and only graphs with that same ``dst`` contribute.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_fundamental_slide"></a>

#### transition\_fundamental\_slide

```python
@classmethod
def transition_fundamental_slide(cls, key)
```

Expand a word in `FundamentalSlideBasis` as the transpose of the polynomial side: the
coefficient of ``candidate`` is the coefficient of the monomial ``x^key`` in the fundamental
slide polynomial of ``candidate``, over all weak compositions of the same length and size.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
@cache
def transition_grothendieck(cls, key)
```

Transition a word key (composition) to the Grothendieck basis.

Coefficient of ``G_w`` is the number of WC graphs of permutation ``w``
and weight ``key``, multiplied by ``beta^(|key|-inv(w))``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition"></a>

#### transition

```python
@classmethod
@cache
def transition(cls, other_basis)
```

Key -> ``{key: coeff}`` function into ``other_basis``; dispatches to the ``transition_*``
method for each directly supported basis and otherwise goes through `SchubertBasis`.

