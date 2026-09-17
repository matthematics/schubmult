<a id="schubmult.rings.free_algebra.composition_schubert_basis"></a>

# schubmult.rings.free\_algebra.composition\_schubert\_basis

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis"></a>

## CompositionSchubertBasis Objects

```python
class CompositionSchubertBasis(FreeAlgebraBasis)
```

Schubert basis indexed by padded trimcode compositions.

A key is a composition ``c`` and corresponds to Schubert key
``(uncode(c), len(c))``.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list (composition).

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.as_schubert_key"></a>

#### as\_schubert\_key

```python
@classmethod
def as_schubert_key(cls, key)
```

Convert a composition key to a Schubert key ``(Permutation, length)``.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, key)
```

Normalize a key to a composition tuple.

Accepts either a Schubert key ``(Permutation, int)`` or a raw tuple.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Return the composition key for the given RC graph.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.inject"></a>

#### inject

```python
@classmethod
def inject(cls, key1, i, key2, coeff=S.One)
```

Inject *key2* into *key1* at position *i* using Schubert multiplication.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Multiply two composition keys by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.coproduct"></a>

#### coproduct

```python
@classmethod
def coproduct(cls, key)
```

Compute the coproduct by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.bcoproduct"></a>

#### bcoproduct

```python
@classmethod
def bcoproduct(cls, key)
```

Compute the bar-coproduct by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.internal_product"></a>

#### internal\_product

```python
@classmethod
def internal_product(cls, key1, key2, coeff=S.One)
```

Compute the internal product by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.skew_element"></a>

#### skew\_element

```python
@classmethod
def skew_element(cls, w, u, n)
```

Compute the skew element by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual basis (delegates to SchubertBasis).

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from CompositionSchubertBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``CompSchub``-labelled display object for the composition key *k*.

