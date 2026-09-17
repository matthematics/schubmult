<a id="schubmult.rings.polynomial_algebra.monomial_basis"></a>

# schubmult.rings.polynomial\_algebra.monomial\_basis

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis"></a>

## MonomialBasis Objects

```python
class MonomialBasis(PolynomialBasis)
```

Standard monomial basis for the polynomial algebra.

Keys are tuples of nonnegative integers representing exponent vectors.
This is the fundamental basis through which other bases transition
by default, and is dual to the :class:`WordBasis` of the free algebra.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the deconcatenation coproduct on a monomial key.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two monomial keys by component-wise addition of exponents.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.expand_monom"></a>

#### expand\_monom

```python
def expand_monom(monom)
```

Convert an exponent tuple to a monomial expression in the generating set.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a dict of monomial keys into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_slide"></a>

#### transition\_slide

```python
def transition_slide(dct, other_basis)
```

Transition a monomial dict to a slide-type basis via triangular inversion.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_slide_monom"></a>

#### transition\_slide\_monom

```python
def transition_slide_monom(other_basis, monom, coeff=S.One)
```

Express a single monomial in a slide-type basis by dominance-order inversion.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_anti_schubert"></a>

#### transition\_anti\_schubert

```python
def transition_anti_schubert(dct, other_basis)
```

Transition monomials to the anti-Schubert basis by reversing exponents.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct)
```

Transition monomials to the Schubert basis by grouping by length.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_double_forest"></a>

#### transition\_double\_forest

```python
def transition_double_forest(dct, other_basis)
```

Transition monomials to DoubleForestPolyBasis via ForestPolyBasis.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`WordBasis`).

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from monomial basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.from_expr"></a>

#### from\_expr

```python
def from_expr(expr, length=None)
```

Parse a symbolic expression into monomial-basis coefficient dict.

