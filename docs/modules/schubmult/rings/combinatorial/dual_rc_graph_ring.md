<a id="schubmult.rings.combinatorial.dual_rc_graph_ring"></a>

# schubmult.rings.combinatorial.dual\_rc\_graph\_ring

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRingElement"></a>

## DualRCGraphRingElement Objects

```python
class DualRCGraphRingElement(SchubertMonomialRingElement)
```

DualRCGraphRing elements are linear combinations of RCGraph basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
to do this directly.

The product * is well defined for any pair of RC graphs and is the dual product.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRingElement.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(perm)
```

Apply divided difference operator for `perm` to self.
Linear extension of RCGraph.divdiff_perm.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRingElement.divdiff"></a>

#### divdiff

```python
def divdiff(*seq)
```

Sequential divided difference operators.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRing"></a>

## DualRCGraphRing Objects

```python
class DualRCGraphRing(SchubertMonomialRing)
```

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRing.schub"></a>

#### schub

```python
def schub(perm, n=None)
```

Return the DualRCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

