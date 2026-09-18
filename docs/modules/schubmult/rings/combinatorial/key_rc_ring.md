<a id="schubmult.rings.combinatorial.key_rc_ring"></a>

# schubmult.rings.combinatorial.key\_rc\_ring

`KeyRCGraphRing`: `RCGraphRing` quotient modeling key polynomials (Demazure characters).

RC graphs are snapped to a canonical representative of their ``extremal_weight`` class
(matching Edelman-Greene recording tableaux), and ``to_free_algebra_element`` lands in
the dual key basis indexed by ``extremal_weight``.

<a id="schubmult.rings.combinatorial.key_rc_ring.KeyRCGraphRingElement"></a>

## KeyRCGraphRingElement Objects

```python
class KeyRCGraphRingElement(RCGraphRingElement)
```

Element of `KeyRCGraphRing`; converts to the free-algebra key basis via ``extremal_weight``.

<a id="schubmult.rings.combinatorial.key_rc_ring.KeyRCGraphRingElement.to_free_algebra_element"></a>

#### to\_free\_algebra\_element

```python
def to_free_algebra_element(basis=None)
```

Map each RC graph to the dual key basis element indexed by its ``extremal_weight``.

<a id="schubmult.rings.combinatorial.key_rc_ring.KeyRCGraphRing"></a>

## KeyRCGraphRing Objects

```python
class KeyRCGraphRing(RCGraphRing)
```

`RCGraphRing` with products snapped to canonical key-class representatives; see the module docstring.

