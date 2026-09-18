<a id="schubmult.rings.combinatorial.wc_graph_ring"></a>

# schubmult.rings.combinatorial.wc\_graph\_ring

`WCGraphRing`: the ring whose basis elements are `WCGraph`s (the K-theoretic / Grothendieck
analogue of `RCGraphRing`). ``to_free_algebra_element`` lands in the free-algebra Grothendieck
basis by default.

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRingElement"></a>

## WCGraphRingElement Objects

```python
class WCGraphRingElement(SchubertMonomialRingElement)
```

WCGraphRing elements are linear combinations of WCGraph basis elements.

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRing"></a>

## WCGraphRing Objects

```python
class WCGraphRing(SchubertMonomialRing)
```

The ring of `WCGraph`s; see the module docstring.

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRing.groth"></a>

#### groth

```python
def groth(perm, n=None)
```

Return the WCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

