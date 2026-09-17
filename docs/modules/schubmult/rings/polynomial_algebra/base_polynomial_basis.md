<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis"></a>

# schubmult.rings.polynomial\_algebra.base\_polynomial\_basis

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis"></a>

## PolynomialBasis Objects

```python
class PolynomialBasis(ABC)
```

Abstract base class for polynomial algebra bases.

Subclasses define how keys are represented, how to transition between
bases, and how to expand elements into explicit polynomials. Default
implementations delegate through the :class:`MonomialBasis`.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.is_key"></a>

#### is\_key

```python
@abstractmethod
def is_key(x)
```

Return True if *x* is a valid key for this basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.as_key"></a>

#### as\_key

```python
@abstractmethod
def as_key(x)
```

Normalize *x* into a canonical key for this basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.attach_key"></a>

#### attach\_key

```python
def attach_key(dct)
```

Normalize all keys in *dct* via :meth:`as_key`.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.transition"></a>

#### transition

```python
@abstractmethod
def transition(other_basis)
```

Return a function mapping dicts of this basis to dicts in *other_basis*.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.from_expr"></a>

#### from\_expr

```python
def from_expr(expr, length=None)
```

Parse a symbolic expression into this basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.printing_term"></a>

#### printing\_term

```python
@abstractmethod
def printing_term(k)
```

Return the display symbol for key *k*.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.compose_transition"></a>

#### compose\_transition

```python
@staticmethod
def compose_transition(tkeyfunc, output)
```

Apply a transition function to a dict of basis elements.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.change_tensor_basis"></a>

#### change\_tensor\_basis

```python
@classmethod
def change_tensor_basis(cls, tensor_elem, basis1, basis2)
```

Change the bases of both factors of a tensor element.

**Arguments**:

- `tensor_elem` - An element of a tensor product ring.
- `basis1` - Target basis for the left factor.
- `basis2` - Target basis for the right factor.
  

**Returns**:

  The tensor element re-expressed in the new bases.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a basis dict into an explicit polynomial expression.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the coproduct of *key* by delegating through the monomial basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two keys by transitioning to the monomial basis and back.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class.

