<a id="schubmult.rings.combinatorial.crystal_graph_ring"></a>

# schubmult.rings.combinatorial.crystal\_graph\_ring

`CrystalGraphRing`: a ring whose basis elements are crystal-graph objects, with the crystal
operators extended linearly to ring elements. Base class for the RC graph / WC graph /
factor-algebra rings in this package.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRing"></a>

## CrystalGraphRing Objects

```python
class CrystalGraphRing(BaseRing)
```

Ring whose basis elements are CrystalGraph-like objects.

We deliberately do not special-case tensor objects here: CrystalGraphTensor
implements the same CrystalGraph API and will be handled by polymorphism.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRing.dtype"></a>

#### dtype

```python
def dtype()
```

A fresh empty element bound to this ring.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement"></a>

## CrystalGraphRingElement Objects

```python
class CrystalGraphRingElement(BaseRingElement, CrystalGraph)
```

Element of the CrystalGraphRing.

Keys are arbitrary objects that implement the CrystalGraph API (including
CrystalGraphTensor). All crystal operators / statistics are lifted linearly
by delegating to the underlying key's methods.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.phi"></a>

#### phi

```python
def phi(index: int) -> int
```

Maximum of ``phi(index)`` over the basis keys.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.epsilon"></a>

#### epsilon

```python
def epsilon(index: int) -> int
```

Maximum of ``epsilon(index)`` over the basis keys.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index: int)
```

Linearized raising operator: delegate to each key's raising_operator
and collect results in the ring.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index: int)
```

Linearized lowering operator: delegate to each key's lowering_operator.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.crystal_length"></a>

#### crystal\_length

```python
def crystal_length() -> int
```

Maximum of ``crystal_length()`` over the basis keys.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight() -> Tuple["CrystalGraphRingElement", Tuple[int, ...]]
```

Apply linearized raising operators until none changes the element; returns ``(element, raise_seq)``.

