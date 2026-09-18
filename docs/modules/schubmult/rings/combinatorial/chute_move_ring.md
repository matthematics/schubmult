<a id="schubmult.rings.combinatorial.chute_move_ring"></a>

# schubmult.rings.combinatorial.chute\_move\_ring

`ChuteMoveRing`: a `SchubertMonomialRing` whose basis elements are `ChuteMoveElement`s
(RC graphs marked with a set of simultaneous chute-move rows).

<a id="schubmult.rings.combinatorial.chute_move_ring.ChuteMoveRingElement"></a>

## ChuteMoveRingElement Objects

```python
class ChuteMoveRingElement(SchubertMonomialRingElement)
```

ChuteMoveRing elements are linear combinations of ChuteMoveElement basis elements.

<a id="schubmult.rings.combinatorial.chute_move_ring.ChuteMoveRing"></a>

## ChuteMoveRing Objects

```python
class ChuteMoveRing(SchubertMonomialRing)
```

The ring of `ChuteMoveElement`s; products use `ChuteMoveElement.product`.

