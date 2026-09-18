<a id="schubmult.rings.free_algebra.separated_descents_basis"></a>

# schubmult.rings.free\_algebra.separated\_descents\_basis

`SeparatedDescentsBasis(k)`: level-``k`` refinements of `SchubertBasis` in which a key
``(u, v, numvars)`` splits the Schubert index into a factor ``u`` and a factor ``v`` whose
descents are separated at ``k`` (the basis dual to the separated-descents factorization of
`schubmult.rings.schubert.separated_descents`).

`SeparatedDescentsBasis` is a factory producing a ``_SeparatedDescentsBasis`` subclass with
class attribute ``k``; products go through `SchubertBasis`, and `SchubertBasis` expands into
this basis via `SchubertBasis.transition_separated_descents`.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis"></a>

## \_SeparatedDescentsBasis Objects

```python
class _SeparatedDescentsBasis(FreeAlgebraBasis)
```

Separated descents basis of the free algebra (parameterized by level *k*).

Keys are ``(Permutation, Permutation, int)`` triples representing a
factorization into descents above and below a cutoff level *k*.
Instances are created by the :func:`SeparatedDescentsBasis` factory.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a valid separated descents key.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(Permutation, Permutation, int)`` key.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Multiply two separated descents keys via the Schubert basis.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, perm0, perm1, numvars)
```

Transition a separated descents key to the Schubert basis.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, perm0, perm1, n)
```

Transition a separated descents key to the word basis via the Schubert basis.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from this separated descents basis to *other_basis*.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``SepDesc<k>``-prefixed symbol for key *k*.

<a id="schubmult.rings.free_algebra.separated_descents_basis.SeparatedDescentsBasis"></a>

#### SeparatedDescentsBasis

```python
def SeparatedDescentsBasis(k)
```

Factory that creates a separated descents basis class for level *k*.

