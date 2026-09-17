<a id="schubmult.rings.combinatorial.chute_move_ring"></a>

# schubmult.rings.combinatorial.chute\_move\_ring

<a id="schubmult.rings.combinatorial.chute_move_ring.ChuteMoveRingElement"></a>

## ChuteMoveRingElement Objects

```python
class ChuteMoveRingElement(SchubertMonomialRingElement)
```

ChuteMoveRing elements are linear combinations of ChuteMoveElement basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
to do this directly.

The product * is well defined for any pair of RC graphs and is the dual product.

