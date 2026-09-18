<a id="schubmult.rings.free_algebra.schubert_schur_basis"></a>

# schubmult.rings.free\_algebra.schubert\_schur\_basis

`SchubertSchurBasis`: free-algebra basis dual to products ``s_lambda(x_1..x_n) * S_perm``
of a Schur polynomial in the first ``n`` variables with a Schubert polynomial. Keys are
``(partition, perm, numvars)``.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis"></a>

## SchubertSchurBasis Objects

```python
class SchubertSchurBasis(FreeAlgebraBasis)
```

Schubert-Schur basis of the free algebra.

Keys are ``(partition_tuple, Permutation)`` pairs encoding products of
a Grassmannian Schubert polynomial with a Schur polynomial.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a ``(list/tuple, Permutation/list/tuple)`` pair.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(tuple, Permutation)`` key.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Compute the coproduct of a Schubert-Schur key via the Schubert basis.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, lambd, perm)
```

Transition a Schubert-Schur key ``(lambda, perm)`` to the Schubert basis.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, lambd, perm)
```

Transition a Schubert-Schur key to the word basis via the Schubert basis.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from SchubertSchurBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``SS``-prefixed symbol for key *k*.

