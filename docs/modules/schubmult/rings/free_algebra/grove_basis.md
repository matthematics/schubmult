<a id="schubmult.rings.free_algebra.grove_basis"></a>

# schubmult.rings.free\_algebra.grove\_basis

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis"></a>

## GroveBasis Objects

```python
class GroveBasis(FreeAlgebraBasis)
```

Grove basis of the free algebra.

Keys are tuples representing indexed-grove weight vectors.
Transitions to the Grothendieck basis use RC graph enumeration
filtered by grove weight.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``GroveDual``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the GrovePolyBasis as the dual of GroveBasis.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
def transition_grothendieck(cls, key)
```

Transition a grove key to the Grothendieck basis via WC graph enumeration.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from GroveBasis to *other_basis*.

