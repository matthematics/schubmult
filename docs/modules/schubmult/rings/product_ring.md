<a id="schubmult.rings.product_ring"></a>

# schubmult.rings.product\_ring

`ProductRing`: array-backed componentwise product of Schubert-family rings (draft).

An element holds a numpy object array ``_arr`` with one ring element per factor; arithmetic is
componentwise on the array. This is an older sketch of the same idea as
`schubmult.rings.direct_product_ring.DirectProductRing`, which is the supported implementation;
`ProductRing` is not exported from the package.

<a id="schubmult.rings.product_ring.ProductRing"></a>

## ProductRing Objects

```python
class ProductRing(BaseSchubertRing)
```

Componentwise product of rings backed by a numpy array of factor elements. See the module docstring.

<a id="schubmult.rings.product_ring.ProductRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(*rings)
```

Flatten nested product rings and pool the factors' generators and coefficient generators.

<a id="schubmult.rings.product_ring.ProductRing.rings"></a>

#### rings

```python
@property
def rings()
```

The (flattened) tuple of factor rings.

<a id="schubmult.rings.product_ring.ProductRing.new"></a>

#### new

```python
def new(x)
```

Wrap a sequence of factor elements as an element of this ring.

<a id="schubmult.rings.product_ring.ProductRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(*x)
```

Build an element from one input per factor (or a single sequence/element), coercing each
input in its factor ring.

<a id="schubmult.rings.product_ring.ProductBasisElement"></a>

## ProductBasisElement Objects

```python
class ProductBasisElement(PrintingTerm)
```

Printing term for a `ProductRing` element; renders the factors joined by ``#``.

<a id="schubmult.rings.product_ring.ProductRingElement"></a>

## ProductRingElement Objects

```python
class ProductRingElement(BaseSchubertElement)
```

Element of a `ProductRing`; arithmetic operators act componentwise on ``_arr``.

