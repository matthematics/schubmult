<a id="schubmult.rings.free_algebra.elementary_basis"></a>

# schubmult.rings.free\_algebra.elementary\_basis

`ElementaryBasis`: free-algebra basis indexed by ``(composition, numvars)``, dual to products of
elementary symmetric polynomials ``e_{a_1}(x_1) e_{a_2}(x_1, x_2) ... e_{a_{n-1}}(x_1..x_{n-1})`` times a
symmetric tail of ``e_k(x_1..x_n)`` factors (`ElemSymPolyBasis`).

Transitions to and from `SchubertBasis` go through the finite ``(numvars, degree)`` block: expanding
every elementary product of that degree in the Schubert basis (Pieri rule) gives the Schubert -> Elem
matrix, and its inverse is Elem -> Schubert.

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

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.canonical_key"></a>

#### canonical\_key

```python
@classmethod
def canonical_key(cls, tup, numvars)
```

Canonicalize ``(tup, numvars)``: the flag part ``tup[:numvars-1]`` is kept as is; the
symmetric tail ``tup[numvars-1:]`` (indices of ``e_k(x_1..x_numvars)`` factors) is sorted
with zeros dropped, or ``(0,)`` if empty.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.degree_keys"></a>

#### degree\_keys

```python
@staticmethod
def degree_keys(numvars, degree)
```

All canonical keys of the given degree in ``numvars`` variables: flag part ``a_i <= i``,
tail a partition with parts in ``1..numvars``. There are as many as monomials of that degree.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.schubert_block"></a>

#### schubert\_block

```python
@classmethod
@cache
def schubert_block(cls, numvars, degree)
```

``(keys, perms, to_schubert, to_elementary)`` for one ``(numvars, degree)`` block.

``perms`` are the permutations whose Schubert polynomial lies in ``Z[x_1..x_numvars]`` with
that degree (Lehmer codes of length ``numvars``). ``to_schubert[key][perm]`` is the
coefficient of ``S_perm`` in ``E_key``; ``to_elementary[perm][key]`` is the inverse matrix,
i.e. the coefficient of ``Elem(key)`` in ``Schub(perm)``.

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

Transition an elementary key to the Schubert basis (row of the inverse block matrix).

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.transition_from_schubert"></a>

#### transition\_from\_schubert

```python
@classmethod
def transition_from_schubert(cls, perm, numvars)
```

Expand ``Schub(perm, numvars)`` in this basis: the Schubert coefficients of each ``E_key``.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``Elem``-labelled display object for key *k*.

