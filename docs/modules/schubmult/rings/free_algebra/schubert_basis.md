<a id="schubmult.rings.free_algebra.schubert_basis"></a>

# schubmult.rings.free\_algebra.schubert\_basis

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis"></a>

## SchubertBasis Objects

```python
class SchubertBasis(FreeAlgebraBasis)
```

Schubert basis of the free algebra.

Keys are ``(Permutation, int)`` pairs where the integer records the
number of variables.  Products use the separated-descents Schubert
ring, and transitions to the word basis go through elementary
symmetric function decompositions.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a valid Schubert basis key.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Return the Schubert key ``(perm, length)`` for the given RC graph.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(Permutation, numvars)`` key.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply two Schubert basis keys via the separated-descents ring.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.skew_element"></a>

#### skew\_element

```python
@classmethod
def skew_element(cls, w, u, n)
```

Compute the skew Schubert element S_w / S_u truncated to *n* variables.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Compute the coproduct of a Schubert key in the tensor ring.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
@cache
def transition_grothendieck(cls, perm, numvars)
```

Transition a Schubert key to the Grothendieck basis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_schubert_schur"></a>

#### transition\_schubert\_schur

```python
@classmethod
def transition_schubert_schur(cls, *x)
```

Transition a Schubert key to the Schubert-Schur basis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_schur_elementary"></a>

#### transition\_schur\_elementary

```python
@classmethod
def transition_schur_elementary(cls, *x)
```

Transition a Schubert key to the Schur-Elementary basis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_elementary"></a>

#### transition\_elementary

```python
@classmethod
def transition_elementary(cls, perm, numvars)
```

Transition a Schubert key to the elementary basis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_separated_descents"></a>

#### transition\_separated\_descents

```python
@classmethod
def transition_separated_descents(cls, k, *x)
```

Transition a Schubert key to the separated descents basis of level *k*.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_jbasis"></a>

#### transition\_jbasis

```python
@classmethod
def transition_jbasis(cls, perm, n)
```

Transition a Schubert key ``(perm, n)`` to the J basis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the SchubertPolyBasis as the dual of SchubertBasis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition"></a>

#### transition

```python
@classmethod
@cache
def transition(cls, other_basis)
```

Return a transition function from SchubertBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.old_transition_word"></a>

#### old\_transition\_word

```python
@classmethod
@cache
def old_transition_word(cls, perm, numvars)
```

Transition to WordBasis via SEM basis (legacy implementation).

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
@cache
def transition_word(cls, perm, numvars)
```

Transition a Schubert key to the word basis via SEM factorization.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a Xi-notation display object for key *k*.

