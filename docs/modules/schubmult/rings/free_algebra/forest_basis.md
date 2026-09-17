<a id="schubmult.rings.free_algebra.forest_basis"></a>

# schubmult.rings.free\_algebra.forest\_basis

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis"></a>

## ForestBasis Objects

```python
class ForestBasis(FreeAlgebraBasis)
```

Forest basis of the free algebra.

Keys are tuples representing indexed-forest weight vectors.
Transitions to the Schubert basis use RC graph enumeration
filtered by forest weight.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``ForestDual``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the ForestPolyBasis as the dual of ForestBasis.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a forest key to the Schubert basis via RC graph enumeration.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from ForestBasis to *other_basis*.

