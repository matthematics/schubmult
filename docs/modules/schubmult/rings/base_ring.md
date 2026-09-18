<a id="schubmult.rings.base_ring"></a>

# schubmult.rings.base\_ring

Shared machinery for every ring in `schubmult.rings`.

`BaseRingElement` is a ``dict`` mapping basis keys (permutations, RC graphs, tuples,
...) to coefficients, wired into sympy's printing and arithmetic protocols so that
elements can be added, multiplied, and displayed. `BaseRing` provides the
corresponding ring-level operations (``add``/``sub``/``mul``, coercion via
``domain_new``, construction via ``from_dict``/``from_expr``) and declares the hooks
concrete rings must implement (``new``, ``printing_term``, ``mul_expr``, ...).
A ring's element type is created dynamically as ``self.dtype`` with ``ring`` bound.

<a id="schubmult.rings.base_ring.BaseRingElement"></a>

## BaseRingElement Objects

```python
class BaseRingElement(DomainElement, DefaultPrinting, dict)
```

A ring element: ``{basis_key: coefficient}`` with sympy-compatible arithmetic and printing.

<a id="schubmult.rings.base_ring.BaseRingElement.is_zero"></a>

#### is\_zero

```python
@property
def is_zero()
```

Whether every coefficient is exactly zero.

<a id="schubmult.rings.base_ring.BaseRingElement.parent"></a>

#### parent

```python
def parent()
```

The ring this element belongs to (sympy domain protocol).

<a id="schubmult.rings.base_ring.BaseRingElement.has_free"></a>

#### has\_free

```python
def has_free(*args)
```

Whether any of the given symbols appears in ``free_symbols``.

<a id="schubmult.rings.base_ring.BaseRingElement.apply_to_keys"></a>

#### apply\_to\_keys

```python
def apply_to_keys(func)
```

Map each basis key through ``func`` (dropping keys where it returns ``None``), keeping coefficients.

<a id="schubmult.rings.base_ring.BaseRingElement.as_terms"></a>

#### as\_terms

```python
def as_terms()
```

Terms ``coeff * basis_symbol`` in dict order (sympy printing hook).

<a id="schubmult.rings.base_ring.BaseRingElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms sorted by basis key (sympy printing hook).

<a id="schubmult.rings.base_ring.BaseRingElement.coproduct"></a>

#### coproduct

```python
def coproduct()
```

Coproduct into the tensor square ring, via ``ring.coproduct_on_basis``.

<a id="schubmult.rings.base_ring.BaseRingElement.as_coefficients_dict"></a>

#### as\_coefficients\_dict

```python
def as_coefficients_dict()
```

``{basis_symbol: coeff}`` mapping display symbols to coefficients.

<a id="schubmult.rings.base_ring.BaseRingElement.expand"></a>

#### expand

```python
def expand(deep=True, *args, **kwargs)
```

With ``deep=True`` expand to an explicit polynomial (``as_polynomial``); with ``deep=False`` only
expand each coefficient, keeping the basis.

<a id="schubmult.rings.base_ring.BaseRingElement.as_expr"></a>

#### as\_expr

```python
def as_expr()
```

Sum of the ``as_terms()`` as a sympy ``Add``.

<a id="schubmult.rings.base_ring.BaseRingElement.as_polynomial"></a>

#### as\_polynomial

```python
def as_polynomial()
```

Hook: expand this element to an explicit polynomial expression.

<a id="schubmult.rings.base_ring.BaseRingElement.almosteq"></a>

#### almosteq

```python
def almosteq(other)
```

Equality up to coefficient expansion.

<a id="schubmult.rings.base_ring.BaseRingElement.__matmul__"></a>

#### \_\_matmul\_\_

```python
def __matmul__(other)
```

Tensor product ``self (x) other`` in the `TensorRing` of the two rings.

<a id="schubmult.rings.base_ring.BaseRing"></a>

## BaseRing Objects

```python
class BaseRing(Ring, CompositeDomain)
```

Abstract base ring over a sympy coefficient domain (default ``EXRAW``); see the module docstring.

<a id="schubmult.rings.base_ring.BaseRing.__matmul__"></a>

#### \_\_matmul\_\_

```python
def __matmul__(other)
```

The `TensorRing` ``self (x) other``.

<a id="schubmult.rings.base_ring.BaseRing.to_sympy"></a>

#### to\_sympy

```python
def to_sympy(elem)
```

Convert an element to a sympy expression (``as_expr``).

<a id="schubmult.rings.base_ring.BaseRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(domain=None)
```

**Arguments**:

- `domain` - Coefficient domain; defaults to sympy's ``EXRAW`` (arbitrary expressions).

<a id="schubmult.rings.base_ring.BaseRing.add"></a>

#### add

```python
def add(elem, other)
```

Coefficient-wise sum, dropping zeros.

<a id="schubmult.rings.base_ring.BaseRing.sub"></a>

#### sub

```python
def sub(elem, other)
```

Coefficient-wise difference, dropping zeros.

<a id="schubmult.rings.base_ring.BaseRing.neg"></a>

#### neg

```python
def neg(elem)
```

Negate every coefficient.

<a id="schubmult.rings.base_ring.BaseRing.rmul"></a>

#### rmul

```python
def rmul(elem, other)
```

Right-multiply by a scalar (via ``domain_new``), falling back to ``mul_expr``.

<a id="schubmult.rings.base_ring.BaseRing.mul"></a>

#### mul

```python
def mul(elem, other)
```

Multiply by a scalar (via ``domain_new``), falling back to ``mul_expr``.

<a id="schubmult.rings.base_ring.BaseRing.from_sympy"></a>

#### from\_sympy

```python
def from_sympy(expr)
```

Alias for ``from_expr`` (sympy domain protocol).

<a id="schubmult.rings.base_ring.BaseRing.new"></a>

#### new

```python
def new(x)
```

Hook: build an element from ``x``.

<a id="schubmult.rings.base_ring.BaseRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Hook: the sympy symbol displayed for basis key ``k``.

<a id="schubmult.rings.base_ring.BaseRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
def coproduct_on_basis(k)
```

Hook: coproduct of basis key ``k`` in the tensor square ring.

<a id="schubmult.rings.base_ring.BaseRing.one"></a>

#### one

```python
@property
def one()
```

The multiplicative identity: coefficient 1 on ``zero_monom``.

<a id="schubmult.rings.base_ring.BaseRing.is_elem_mul_type"></a>

#### is\_elem\_mul\_type

```python
def is_elem_mul_type(elem)
```

Hook: whether ``elem`` should use the ``elem_mul`` fast path.

<a id="schubmult.rings.base_ring.BaseRing.elem_mul"></a>

#### elem\_mul

```python
def elem_mul(ring_elem, elem)
```

Hook: elementary-symmetric fast-path multiplication.

<a id="schubmult.rings.base_ring.BaseRing.from_dict"></a>

#### from\_dict

```python
def from_dict(element, orig_domain=None)
```

Build an element from ``{key: coeff}``, coercing each coefficient via ``domain_new`` and dropping zeros.

<a id="schubmult.rings.base_ring.BaseRing.from_dict_unchecked"></a>

#### from\_dict\_unchecked

```python
def from_dict_unchecked(element)
```

from_dict for coefficients already known to lie in the domain (drops structural zeros only).

<a id="schubmult.rings.base_ring.BaseRing.zero"></a>

#### zero

```python
@property
def zero()
```

The empty element.

<a id="schubmult.rings.base_ring.BaseRing.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce ``element`` into the coefficient domain (``sympify``), refusing ring/domain elements.

<a id="schubmult.rings.base_ring.BaseRing.from_expr"></a>

#### from\_expr

```python
def from_expr(x)
```

Build an element from an expression by multiplying the identity by it.

<a id="schubmult.rings.base_ring.BaseRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Hook: multiply ``elem`` by a symbolic expression ``x``.

