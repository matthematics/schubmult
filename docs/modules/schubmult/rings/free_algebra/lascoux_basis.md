<a id="schubmult.rings.free_algebra.lascoux_basis"></a>

# schubmult.rings.free\_algebra.lascoux\_basis

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis"></a>

## LascouxBasis Objects

```python
class LascouxBasis(FreeAlgebraBasis)
```

Lascoux polynomial (Demazure character) basis of the free algebra.

Lascouxs are weak composition tuples. Transitions to the Schubert basis
use RC graph enumeration filtered by extremal weight.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the LascouxPolyBasis as the dual of LascouxBasis.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``Lascoux``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
def transition_grothendieck(cls, key)
```

Transition a Lascoux composition to the Grothendieck basis via WC graphs.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from LascouxBasis to *other_basis*.

