<a id="schubmult.rings.polynomial_algebra._core"></a>

# schubmult.rings.polynomial\_algebra.\_core

`PolynomialAlgebra`: the polynomial ring ``Z[x_1, x_2, ...]`` with a pluggable basis.

The ring itself is basis-agnostic; a `PolynomialBasis` instance supplies the key
type, the product rule, the coproduct, and the transitions to/from the monomial
basis. Elements of rings with different bases are interconverted via
``change_basis``. ``PA`` is the standard monomial-basis instance in ``x``; the
pre-built instances for other bases (``Schub``, ``Key``, ``FSlide``, ...) live in
the package ``__init__``.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement"></a>

## PolynomialAlgebraElement Objects

```python
class PolynomialAlgebraElement(BaseRingElement)
```

Element of a polynomial algebra, stored as a dict mapping basis keys to coefficients.

Keys are exponent tuples (in the monomial basis) or basis-specific keys
depending on the parent ring's basis. Supports arithmetic, basis changes,
and duality pairing with free algebra elements.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.as_coefficients_dict"></a>

#### as\_coefficients\_dict

```python
def as_coefficients_dict()
```

Return a dict mapping printing terms to sympified coefficients.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.branch"></a>

#### branch

```python
def branch(index)
```

Split the variables at ``index``: ``x_1..x_index`` on the left tensor factor, the rest on the
right, returned in the tensor square of this ring's basis.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.coproduct"></a>

#### coproduct

```python
def coproduct()
```

Sum of ``branch(index)`` over every split point (the full variable-splitting coproduct).

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.change_basis"></a>

#### change\_basis

```python
def change_basis(other_basis: type)
```

Convert this element to another polynomial basis.

**Arguments**:

- `other_basis` - A basis class, basis instance, or callable returning a basis.
  

**Returns**:

  A new PolynomialAlgebraElement in the target basis's ring.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.expand"></a>

#### expand

```python
def expand()
```

Expand this element into an explicit polynomial expression.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.apply_dual_element"></a>

#### apply\_dual\_element

```python
def apply_dual_element(dual_elem)
```

Pair this polynomial element with a dual free algebra element.

Converts *self* to the monomial basis and *dual_elem* to the word
basis, then sums products of matching coefficients.

**Arguments**:

- `dual_elem` - A FreeAlgebraElement to pair with.
  

**Returns**:

  The scalar pairing value.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra"></a>

## PolynomialAlgebra Objects

```python
class PolynomialAlgebra(BaseRing)
```

Polynomial algebra ring with a configurable basis.

The algebra operates on :class:`PolynomialAlgebraElement` instances whose
keys are determined by the chosen basis. Supports multiplication, basis
changes, coproducts, and conversion from symbolic expressions.

**Arguments**:

- `basis` - A basis instance (e.g. ``MonomialBasis(x)``).
- `domain` - Coefficient domain (default ``EXRAW``).

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.__init__"></a>

#### \_\_init\_\_

```python
def __init__(basis, domain=None)
```

Initialize a PolynomialAlgebra with the given basis and coefficient domain.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.genset"></a>

#### genset

```python
@property
def genset()
```

The basis's generating set.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
@cache
def coproduct_on_basis(key)
```

Compute the coproduct of a single basis key in the tensor ring.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.mul"></a>

#### mul

```python
def mul(elem, other)
```

Multiply two elements via the basis product rule.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.new"></a>

#### new

```python
def new(*x)
```

Create a new element from the given key or expression.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.from_expr"></a>

#### from\_expr

```python
def from_expr(x, length=None)
```

Create an element from a symbolic expression.

Parses *x* into monomials, then transitions to this ring's basis.

**Arguments**:

- `x` - A symbolic polynomial expression.
- `length` - Optional fixed number of variables.
  

**Returns**:

  A PolynomialAlgebraElement in this ring.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Return the display symbol for basis key *k*.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.from_dict"></a>

#### from\_dict

```python
def from_dict(element)
```

Construct an element from a dict of ``{key: coefficient}`` pairs.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce a raw value into the coefficient domain.

