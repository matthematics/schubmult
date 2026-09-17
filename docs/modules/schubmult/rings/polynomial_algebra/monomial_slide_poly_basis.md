<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.monomial\_slide\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis"></a>

## MonomialSlidePolyBasis Objects

```python
class MonomialSlidePolyBasis(PolynomialBasis)
```

Monomial slide polynomial basis.

Keys are weak compositions. Monomial slide polynomials refine
key polynomials and are coarser than monomials, using a recursive
construction based on the first nonzero entry.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a monomial slide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from monomial slide basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a monomial slide basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from monomial slide basis to *other_basis*.

