<a id="schubmult.rings.free_algebra.glide_basis"></a>

# schubmult.rings.free\_algebra.glide\_basis

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis"></a>

## GlideBasis Objects

```python
class GlideBasis(FreeAlgebraBasis)
```

Glide basis of the free algebra.

Keys are weak composition tuples. Transitions are computed via
polynomial algebra slide polynomial expansions.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``FS``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
def transition_grothendieck(cls, key)
```

Transition a glide key to the Grothendieck basis.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the GlidePolyBasis as the dual.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from GlideBasis to *other_basis*.

