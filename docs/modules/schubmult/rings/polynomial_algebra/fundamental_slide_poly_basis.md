<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.fundamental\_slide\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.get_descent_composition"></a>

#### get\_descent\_composition

```python
def get_descent_composition(word)
```

Compute the descent composition of a word.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.slide_product"></a>

#### slide\_product

```python
def slide_product(a, b)
```

Compute the structure constants for multiplying two fundamental slide polynomials.

Given weak compositions *a* and *b*, returns a dict mapping result
compositions to their coefficients in the fundamental slide expansion
of the product.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis"></a>

## FundamentalSlidePolyBasis Objects

```python
class FundamentalSlidePolyBasis(PolynomialBasis)
```

Fundamental slide polynomial basis.

Keys are weak compositions. Fundamental slide polynomials provide a
basis that refines Schubert polynomials and coarsens monomials, with
an efficient combinatorial product rule.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a slide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`FundamentalSlideBasis`).

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a slide basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from fundamental slide basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from fundamental slide basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two fundamental slide keys using the slide product rule.

