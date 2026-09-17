<a id="schubmult.rings.free_algebra.schur_elementary_basis"></a>

# schubmult.rings.free\_algebra.schur\_elementary\_basis

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis"></a>

## SchurElementaryBasis Objects

```python
class SchurElementaryBasis(FreeAlgebraBasis)
```

Schur-Elementary basis of the free algebra.

Keys are ``(tuple, tuple)`` pairs encoding products of
a standard elementary monomial (first tuple) and a Schur polynomial (second tuple, partition
of length precisely len(first_tuple) + 1 in increasing order, with zeros at the beginning if
needed).  For example, the key ``((1, 2), (0, 1, 3))`` corresponds to the product of the elementary monomial

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a ``(list/tuple, list/tuple)`` pair.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(tuple, tuple)`` key.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Compute the coproduct of a Schubert-Schur key via the Schubert basis.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, elem_tup, lambd)
```

Transition a Schubert-Schur key ``(lambda, perm)`` to the Schubert basis.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, elem_tup, lambd)
```

Transition a Schubert-Schur key to the word basis via the Schubert basis.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from SchurElementaryBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``SE``-prefixed symbol for key *k*.

