<a id="schubmult.rings.combinatorial.bpd_ring"></a>

# schubmult.rings.combinatorial.bpd\_ring

`BPDRing`: a `SchubertMonomialRing` whose basis elements are bumpless pipe dreams (`BPD`),
with conversion to `RCGraphRing` via ``to_rc_graph_ring_element``.

<a id="schubmult.rings.combinatorial.bpd_ring.BPDRingElement"></a>

## BPDRingElement Objects

```python
class BPDRingElement(SchubertMonomialRingElement)
```

Linear combination of `BPD` basis elements.

<a id="schubmult.rings.combinatorial.bpd_ring.BPDRingElement.to_rc_graph_ring_element"></a>

#### to\_rc\_graph\_ring\_element

```python
def to_rc_graph_ring_element(
        rc_ring: RCGraphRing | None = None) -> RCGraphRingElement
```

Convert each BPD to its RC graph and re-express in an `RCGraphRing`.

<a id="schubmult.rings.combinatorial.bpd_ring.BPDRing"></a>

## BPDRing Objects

```python
class BPDRing(SchubertMonomialRing)
```

The ring of bumpless pipe dreams; products use `BPD.product`.

