<a id="schubmult.rings.tensor_ring"></a>

# schubmult.rings.tensor\_ring

`TensorRing`: tensor products ``R_1 (x) ... (x) R_n`` of `BaseRing` instances.

Built with the ``@`` operator on rings (``Sx @ Sx``) or ``TensorRing(R1, R2, ...)``; nested
tensor rings are flattened. Keys are tuples ``(k_1, ..., k_n)`` of factor keys, multiplication
is factorwise, and the coproduct of a ring lands in ``R @ R``. Elements print as
``a # b``.

<a id="schubmult.rings.tensor_ring.TensorRing"></a>

## TensorRing Objects

```python
class TensorRing(BaseRing)
```

Tensor product of rings; keys are tuples of factor keys. See the module docstring.

<a id="schubmult.rings.tensor_ring.TensorRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
def coproduct_on_basis(k)
```

Coproduct of the basis element ``k = (k_1, ..., k_n)`` into ``self @ self``.

Takes each factor's coproduct and interlaces them into flat keys
``(k_1^L, ..., k_n^L, k_1^R, ..., k_n^R)``.

<a id="schubmult.rings.tensor_ring.TensorRing.from_rc_graph_tensor"></a>

#### from\_rc\_graph\_tensor

```python
def from_rc_graph_tensor(rc_graph_tensor)
```

Pure tensor of the two factor rings' ``from_rc_graph`` images of a pair of RC graphs.

<a id="schubmult.rings.tensor_ring.TensorRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(*rings)
```

Tensor the given rings, flattening any that are themselves tensor rings.

<a id="schubmult.rings.tensor_ring.TensorRing.rings"></a>

#### rings

```python
@property
def rings()
```

The (flattened) tuple of tensor factors.

<a id="schubmult.rings.tensor_ring.TensorRing.rmul"></a>

#### rmul

```python
def rmul(elem1, elem2)
```

Scale every coefficient of ``elem1`` by the scalar ``elem2``.

<a id="schubmult.rings.tensor_ring.TensorRing.mul"></a>

#### mul

```python
def mul(elem1, elem2)
```

Factorwise product: ``(a_1 (x) ... (x) a_n) * (b_1 (x) ... (x) b_n) = (a_1 b_1) (x) ... (x) (a_n b_n)``,
expanding each factor product in its own ring.

<a id="schubmult.rings.tensor_ring.TensorRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

Product of the factor rings' polynomials for the key tuple ``k``.

<a id="schubmult.rings.tensor_ring.TensorRing.from_comp_ring"></a>

#### from\_comp\_ring

```python
def from_comp_ring(t)
```

Embed an element of one factor (or of a sub-tensor of factors) into this ring, filling the
other positions with their ``zero_monom`` (the identity).

<a id="schubmult.rings.tensor_ring.TensorRing.ext_multiply"></a>

#### ext\_multiply

```python
def ext_multiply(elem1, elem2)
```

External (tensor) product ``elem1 (x) elem2``: concatenates keys, flattening tensor-ring inputs.

<a id="schubmult.rings.tensor_ring.TensorRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(x)
```

A key tuple gives the corresponding basis element; anything else is parsed via ``from_expr``.

<a id="schubmult.rings.tensor_ring.TensorBasisElement"></a>

## TensorBasisElement Objects

```python
class TensorBasisElement(PrintingTerm)
```

Printing term for a tensor key; renders as ``a # b`` (str) or a tensor product (pretty/LaTeX).

<a id="schubmult.rings.tensor_ring.TensorRingElement"></a>

## TensorRingElement Objects

```python
class TensorRingElement(BaseRingElement)
```

Element of a `TensorRing`: a dict from key tuples to coefficients.

<a id="schubmult.rings.tensor_ring.TensorRingElement.coproduct"></a>

#### coproduct

```python
def coproduct()
```

Coproduct into ``ring @ ring`` via `TensorRing.coproduct_on_basis`.

<a id="schubmult.rings.tensor_ring.TensorRingElement.expand"></a>

#### expand

```python
def expand(deep=True, *args, **kwargs)
```

Expand to a commutative polynomial by multiplying out the factors' expansions (all factors
are assumed to live in disjoint or commuting variable sets).

