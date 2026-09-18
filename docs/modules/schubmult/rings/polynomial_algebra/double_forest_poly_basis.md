<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.double\_forest\_poly\_basis

`DoubleForestPolyBasis`: the double (two-alphabet) forest polynomial basis of `PolynomialAlgebra`;
see `schubmult.combinatorics.double_forest` for the underlying polynomials.

<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis.DoubleForestPolyBasis"></a>

## DoubleForestPolyBasis Objects

```python
class DoubleForestPolyBasis(PolynomialBasis)
```

Abstract double forest polynomial basis.

Keys are weak compositions indexing double forest basis elements DF[key],
with polynomial coefficients in a second generating set (equivariant vars).

<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis.DoubleForestPolyBasis.basis_polynomial"></a>

#### basis\_polynomial

```python
def basis_polynomial(key)
```

Expand one double-forest basis key as a polynomial in x,t.

<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis.DoubleForestPolyBasis.basis_forest_expansion"></a>

#### basis\_forest\_expansion

```python
def basis_forest_expansion(key, length)
```

Expand one double-forest basis key into ForestPolyBasis in x.

