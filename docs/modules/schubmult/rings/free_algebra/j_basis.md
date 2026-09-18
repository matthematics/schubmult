<a id="schubmult.rings.free_algebra.j_basis"></a>

# schubmult.rings.free\_algebra.j\_basis

`JBasis`: free-algebra basis indexed by compositions with no zeros.

A Schubert key ``(perm, n)`` whose padded code has no zeros is itself a J key; zeros are
handled by the transitions in `SchubertBasis.transition_jbasis` and `WordBasis.transition_jbasis`.

<a id="schubmult.rings.free_algebra.j_basis.JBasis"></a>

## JBasis Objects

```python
class JBasis(FreeAlgebraBasis)
```

J basis of the free algebra.

Keys are tuples of positive integers (no zeros allowed in transitions).
The J basis indexes elements whose Schubert expansion has no zero
entries in the Lehmer code.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.from_perm"></a>

#### from\_perm

```python
@staticmethod
def from_perm(perm, n)
```

Extract a J basis key from *perm* if the first *n* code entries are nonzero.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.coproduct"></a>

#### coproduct

```python
@classmethod
def coproduct(cls, key)
```

Coproduct for JBasis equals the bar-coproduct.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``J``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from JBasis to *other_basis*.

