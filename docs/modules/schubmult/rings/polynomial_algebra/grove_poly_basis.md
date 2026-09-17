<a id="schubmult.rings.polynomial_algebra.grove_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.grove\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis"></a>

## GrovePolyBasis Objects

```python
class GrovePolyBasis(PolynomialBasis)
```

Grove polynomial basis.

Keys are weak compositions encoding indexed forests. Grove polynomials are
the ``beta``-deformed (set-valued) forest polynomials, computed by summing
over set-valued labelings of the corresponding forest structure.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a grove key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a grove basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from grove basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.product"></a>

#### product

```python
@cache
def product(key1, key2, coeff=S.One)
```

Multiply two grove keys using the glide product rule.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from grove basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`GroveBasis`).

