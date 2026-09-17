<a id="schubmult.rings.direct_product_ring"></a>

# schubmult.rings.direct\_product\_ring

<a id="schubmult.rings.direct_product_ring.DirectProductRing"></a>

## DirectProductRing Objects

```python
class DirectProductRing(BaseRing)
```

Direct product of an arbitrary (fixed) number of rings.

An element is stored as a dict mapping keys ``(i, k)`` to coefficients,
where ``i`` is the component index and ``k`` is a basis key from the
*i*-th ring.  Addition and multiplication are componentwise: terms from
different components never interact, and cross-component products are
zero.

Parameters
----------
``*rings`` : BaseRing
    One or more constituent rings.

Examples
--------
>>> D = DirectProductRing(R1, R2, R3)
>>> a = D.from_component(0, some_R1_element)
>>> b = D.from_component(2, some_R3_element)
>>> a + b          # lives in components 0 and 2
>>> a * b          # zero (different components)
>>> D[0]           # R1
>>> D.project(a, 0)  # back to R1

<a id="schubmult.rings.direct_product_ring.DirectProductRing.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(i)
```

Return the *i*-th constituent ring.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.component_one"></a>

#### component\_one

```python
def component_one(i)
```

Return the identity element supported only on component *i*.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.from_component"></a>

#### from\_component

```python
def from_component(i, elem)
```

Lift an element of ``self[i]`` into the direct product.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.project"></a>

#### project

```python
def project(elem, i)
```

Project onto component *i*, returning an element of ``self[i]``.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.mul"></a>

#### mul

```python
def mul(elem1, elem2)
```

Componentwise multiplication.

Only terms in the same component interact; cross-component products
are zero.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(*elems)
```

Construct an element from one element per component.

``D(e0, e1, ..., en)`` lifts each ``ei`` (an element of ``D[i]``)
into the direct product and sums them.

<a id="schubmult.rings.direct_product_ring.DirectProductRingElement"></a>

## DirectProductRingElement Objects

```python
class DirectProductRingElement(BaseRingElement)
```

<a id="schubmult.rings.direct_product_ring.DirectProductRingElement.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(key)
```

Index by component integer or by ``(i, k)`` basis key.

* ``elem[i]`` — project onto component *i* (returns a ``self.ring[i]`` element).
* ``elem[(i, k)]`` — coefficient lookup (standard dict behaviour).

