<a id="schubmult.rings.schubert.grothendieck_ring"></a>

# schubmult.rings.schubert.grothendieck\_ring

Grothendieck polynomial ring (non-equivariant): the ``Gx`` interface.

`GrothendieckRing` is the beta-deformation of `SingleSchubertRing`; basis
elements ``G_w`` are the K-theoretic Schubert classes, with ``beta = 0``
recovering ordinary Schubert polynomials. There is no coefficient alphabet yet
(see `double_grothendieck_ring` for the equivariant version).

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckElement"></a>

## GrothendieckElement Objects

```python
class GrothendieckElement(BaseSchubertElement)
```

Element of a GrothendieckRing, stored as {Permutation: coeff}.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckElement.as_polynomial"></a>

#### as\_polynomial

```python
def as_polynomial()
```

Expand to an explicit polynomial: ``sum coeff * G_w(x)``.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckElement.mult_poly"></a>

#### mult\_poly

```python
def mult_poly(poly)
```

Multiply by an arbitrary polynomial in ``x`` via the Grothendieck Chevalley rule (`mult_poly_groth`).

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing"></a>

## GrothendieckRing Objects

```python
class GrothendieckRing(BaseSchubertRing)
```

Ring of Grothendieck polynomials.

A deformation of the Schubert polynomial ring with parameter beta.
Basis elements G_w satisfy G_u * G_v = sum_w c^w_{u,v}(beta) G_w
where c^w_{u,v}(beta) are polynomials in beta with integer coefficients.
When beta=0, recovers ordinary Schubert polynomials.

Parameters
----------
genset : GeneratingSet
    The generating set (variable alphabet).
beta : sympy/symengine symbol, optional
    The deformation parameter. Defaults to Symbol("β").

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.beta"></a>

#### beta

```python
@property
def beta()
```

The deformation parameter.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

`mult_poly_groth` with this ring's ``beta`` bound.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.from_expr"></a>

#### from\_expr

```python
def from_expr(expr)
```

Convert a polynomial to the Grothendieck basis: expand in Schubert polynomials first, then
change basis Schubert -> Grothendieck.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, expr)
```

Multiply by an expression by first converting it into the Grothendieck basis.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants ``c^w_{u,v}(beta)`` via ``groth_mul_full_with_ring``; only same-ring products supported.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Same as ``cached_product``.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit Grothendieck polynomial ``G_k(x)`` (cached).

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k, prefix="")
```

The ``GrothendieckPoly`` display symbol for basis element ``k``.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list, an element of this ring, or a polynomial expression.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.from_dict"></a>

#### from\_dict

```python
def from_dict(dct)
```

Build an element from ``{Permutation: coeff}``, dropping terms whose coefficient expands to zero.

