<a id="schubmult.rings.free_algebra.key_basis"></a>

# schubmult.rings.free\_algebra.key\_basis

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis"></a>

## KeyBasis Objects

```python
class KeyBasis(FreeAlgebraBasis)
```

Key polynomial (Demazure character) basis of the free algebra.

Keys are weak composition tuples. Transitions to the Schubert basis
use RC graph enumeration filtered by extremal weight.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the KeyPolyBasis as the dual of KeyBasis.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``Key``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a key composition to the Schubert basis via RC graphs.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply two keys by transitioning to WordBasis and back.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from KeyBasis to *other_basis*.

