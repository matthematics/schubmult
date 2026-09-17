<a id="schubmult.rings.free_algebra.word_basis"></a>

# schubmult.rings.free\_algebra.word\_basis

<a id="schubmult.rings.free_algebra.word_basis.WordBasis"></a>

## WordBasis Objects

```python
class WordBasis(FreeAlgebraBasis)
```

Word basis of the free algebra.

Keys are tuples of nonnegative integers representing words. This is the
fundamental basis through which all other bases perform their operations
via basis transitions.

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

Concatenate two words.

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

Compute the additive coproduct of a word.

Decomposes each letter into all (i, key-i) splittings and combines
via a divide-and-conquer tensor product.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.bcoproduct"></a>

#### bcoproduct

```python
@classmethod
@cache
def bcoproduct(cls, key, coeff=S.One)
```

Compute the bar-coproduct of a word.

Like :meth:`coproduct` but drops empty factors (zeros map to
the empty tuple).

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

Compute the internal product of two words via integer matrices (requires SageMath).

Words must not contain zeros. Returns the dict of result
words weighted by *coeff*.

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

