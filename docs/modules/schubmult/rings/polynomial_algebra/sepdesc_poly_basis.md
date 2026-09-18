<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.sepdesc\_poly\_basis

`SepDescPolyBasis`: the separated-descents polynomial basis of `PolynomialAlgebra`, indexed by
``(perm, num_vars)`` pairs (see `schubmult.rings.schubert.separated_descents`).

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis"></a>

## SepDescPolyBasis Objects

```python
class SepDescPolyBasis(PolynomialBasis)
```

Separated descents polynomial basis.

Keys are ``(Permutation, Permutation, k)`` triples. Elements are
products of pairs of Schubert polynomials, parameterized by a
separation level *k*.

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two separated-descents keys by transitioning through Schubert.

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct, other_basis)
```

Transition from separated descents to Schubert basis by multiplying the pair factors.

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from separated descents to *other_basis*.

