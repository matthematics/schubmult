<a id="schubmult.rings.combinatorial.wc_graph_ring"></a>

# schubmult.rings.combinatorial.wc\_graph\_ring

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRingElement"></a>

## WCGraphRingElement Objects

```python
class WCGraphRingElement(SchubertMonomialRingElement)
```

WCGraphRing elements are linear combinations of WCGraph basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
to do this directly.

The product * is well defined for any pair of RC graphs and is the dual product.

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRing"></a>

## WCGraphRing Objects

```python
class WCGraphRing(SchubertMonomialRing)
```

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRing.groth"></a>

#### groth

```python
def groth(perm, n=None)
```

Return the WCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

