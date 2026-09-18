<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.anti\_schubert\_poly\_basis

`AntiSchubertPolyBasis`: the anti-Schubert (``w0``-conjugated Schubert) polynomial basis of `PolynomialAlgebra`.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis"></a>

## AntiSchubertPolyBasis Objects

```python
class AntiSchubertPolyBasis(PolynomialBasis)
```

Anti-Schubert polynomial basis.

Keys are ``(Permutation, length)`` pairs. This basis reverses the
monomial ordering relative to the standard Schubert basis, with
the coproduct correspondingly reversed.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the reversed coproduct of an anti-Schubert key.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two anti-Schubert keys using the underlying Schubert ring.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition_key_key"></a>

#### transition\_key\_key

```python
def transition_key_key(key)
```

Decompose an anti-Schubert polynomial into key polynomials with reversed weights.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition_key"></a>

#### transition\_key

```python
def transition_key(dct)
```

Transition an anti-Schubert dict to the key polynomial basis.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand an anti-Schubert key into reversed monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition_forest"></a>

#### transition\_forest

```python
def transition_forest(dct)
```

Transition an anti-Schubert dict to the forest polynomial basis.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from anti-Schubert basis to *other_basis*.

