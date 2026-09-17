<a id="schubmult.rings.polynomial_algebra.forest_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.forest\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis"></a>

## ForestPolyBasis Objects

```python
class ForestPolyBasis(PolynomialBasis)
```

Forest polynomial basis.

Keys are weak compositions encoding indexed forests. Forest polynomials
are computed by summing over decreasing labelings of the corresponding
forest structure.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a forest key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a forest basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from forest basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.transition_fundamental_slide"></a>

#### transition\_fundamental\_slide

```python
def transition_fundamental_slide(dct)
```

Transition from forest basis to fundamental slide basis.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.to_fundamental_slide"></a>

#### to\_fundamental\_slide

```python
def to_fundamental_slide(key)
```

Express a single forest key in the fundamental slide basis.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`ForestBasis`).

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from forest basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two forest keys by transitioning through the Schubert basis.

