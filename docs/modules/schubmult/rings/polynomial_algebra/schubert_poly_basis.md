<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.schubert\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis"></a>

## SchubertPolyBasis Objects

```python
class SchubertPolyBasis(PolynomialBasis)
```

Schubert polynomial basis.

Keys are ``(Permutation, length)`` pairs. Schubert polynomials form
the canonical basis for the polynomial algebra in Schubert calculus,
dual to the :class:`SchubertBasis` of the free algebra.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the coproduct of a Schubert key by splitting variable sets.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two Schubert keys using the Schubert ring multiplication.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
def transition_grothendieck(dct)
```

Transition a Schubert dict to the Grothendieck polynomial basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_sepdesc"></a>

#### transition\_sepdesc

```python
def transition_sepdesc(dct, other_basis)
```

Transition from Schubert basis to separated descents basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_elementary"></a>

#### transition\_elementary

```python
def transition_elementary(dct, other_basis)
```

Transition from Schubert basis to elementary symmetric basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_key_fundamental_slide"></a>

#### transition\_key\_fundamental\_slide

```python
def transition_key_fundamental_slide(perm, n)
```

Decompose a Schubert polynomial into fundamental slide polynomials via quasi-Yamanouchi RC-graphs.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_fundamental_slide"></a>

#### transition\_fundamental\_slide

```python
def transition_fundamental_slide(dct)
```

Transition a Schubert dict to the fundamental slide basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_key_key"></a>

#### transition\_key\_key

```python
def transition_key_key(key)
```

Decompose a Schubert polynomial into key polynomials via highest-weight RC-graphs.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_key"></a>

#### transition\_key

```python
def transition_key(dct)
```

Transition a Schubert dict to the key polynomial basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a Schubert key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`SchubertBasis`).

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_forest_key"></a>

#### transition\_forest\_key

```python
def transition_forest_key(key)
```

Decompose a Schubert polynomial into forest polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_forest"></a>

#### transition\_forest

```python
def transition_forest(dct)
```

Transition a Schubert dict to the forest polynomial basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from Schubert basis to *other_basis*.

