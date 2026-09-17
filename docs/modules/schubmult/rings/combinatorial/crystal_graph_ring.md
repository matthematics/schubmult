<a id="schubmult.rings.combinatorial.crystal_graph_ring"></a>

# schubmult.rings.combinatorial.crystal\_graph\_ring

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRing"></a>

## CrystalGraphRing Objects

```python
class CrystalGraphRing(BaseRing)
```

Ring whose basis elements are CrystalGraph-like objects.

We deliberately do not special-case tensor objects here: CrystalGraphTensor
implements the same CrystalGraph API and will be handled by polymorphism.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement"></a>

## CrystalGraphRingElement Objects

```python
class CrystalGraphRingElement(BaseRingElement, CrystalGraph)
```

Element of the CrystalGraphRing.

Keys are arbitrary objects that implement the CrystalGraph API (including
CrystalGraphTensor). All crystal operators / statistics are lifted linearly
by delegating to the underlying key's methods.

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

