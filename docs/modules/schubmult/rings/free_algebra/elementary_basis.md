<a id="schubmult.rings.free_algebra.elementary_basis"></a>

# schubmult.rings.free\_algebra.elementary\_basis

`ElementaryBasis`: free-algebra basis indexed by ``(composition, numvars)``, dual to products of
elementary symmetric polynomials ``e_{c_1}(x_1..x_k) e_{c_2}(x_1..x_{k-1}) ...`` in nested
variable sets. `SchubertBasis` expands into it via the monomials of ``S_{perm * w0}``.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis"></a>

## ElementaryBasis Objects

```python
class ElementaryBasis(FreeAlgebraBasis)
```

Elementary symmetric function basis of the free algebra.

Keys are ``(tuple, int)`` pairs where the tuple encodes an elementary
symmetric function composition and the integer is the number of variables.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a ``(tuple/list, int)`` pair.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(tuple, int)`` key.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from ElementaryBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, tup, numvars)
```

Transition an elementary key to the Schubert basis.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``Elem``-labelled display object for key *k*.

