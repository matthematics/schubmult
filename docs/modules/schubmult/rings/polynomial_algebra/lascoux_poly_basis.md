<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.lascoux\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis"></a>

## LascouxPolyBasis Objects

```python
class LascouxPolyBasis(PolynomialBasis)
```

Lascoux polynomial basis.

Keys are weak compositions. Lascoux polynomials provide a
basis that refines Grothendieck polynomials and coarsens monomials, with
an efficient combinatorial product rule.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a glide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`LascouxBasis`).

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a Lascoux basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from Lascoux basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition_glide_key"></a>

#### transition\_glide\_key

```python
def transition_glide_key(key)
```

Transition a Lascoux key to the glide basis.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition_glide"></a>

#### transition\_glide

```python
def transition_glide(dct)
```

Transition from Lascoux basis to glide basis.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from Lascoux basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.product"></a>

#### product

```python
@cache
def product(key1, key2, coeff=S.One)
```

Multiply two Lascoux keys using the Lascoux product rule.

