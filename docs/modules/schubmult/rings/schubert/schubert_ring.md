<a id="schubmult.rings.schubert.schubert_ring"></a>

# schubmult.rings.schubert.schubert\_ring

Ordinary (single) Schubert polynomial ring: the ``Sx`` interface.

`SingleSchubertRing` is a `DoubleSchubertRing` whose coefficient alphabet is
identically zero, so ``S_w(x; 0) = S_w(x)``. Products dispatch to the fast
integer kernel ``schubmult_py`` when both operands are single, and to
``schubmult_double`` when mixed with a genuinely double element.

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing"></a>

## SingleSchubertRing Objects

```python
class SingleSchubertRing(DoubleSchubertRing)
```

The ring of ordinary Schubert polynomials ``S_w(x)``; ``Sx`` is the standard instance.

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants of ``S_u * S_v``: integer ``schubmult_py`` when ``basis2`` is this ring,
else ``schubmult_double`` with ``y = 0``.

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Same as ``cached_product`` (single coefficients are already nonnegative integers).

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.single_variable"></a>

#### single\_variable

```python
def single_variable(elem, varnum)
```

Multiply by ``x_varnum`` (non-equivariant Monk rule).

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list or a polynomial expression.

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.elem_func"></a>

#### elem\_func

```python
@property
def elem_func()
```

`ElemSym` (non-factorial elementary symmetric function).

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.divdiff"></a>

#### divdiff

```python
def divdiff(v, elem)
```

Apply the divided difference ``partial_v``: ``S_u -> S_{u v^{-1}}`` when length-additive, else 0.

