<a id="schubmult.rings.combinatorial.slide_rc_ring"></a>

# schubmult.rings.combinatorial.slide\_rc\_ring

`SlideRCGraphRing`: `RCGraphRing` modeling fundamental slide polynomials.

RC graphs are snapped to a canonical representative of their quasi-Yamanouchi class
(``snap_qy().length_vector``), and ``slide_poly(comp)`` is the sum of RC graphs with
quasi-Yamanouchi weight ``comp``. Products go through `BoundedRCFactorAlgebra`.

<a id="schubmult.rings.combinatorial.slide_rc_ring.SlideRCGraphRingElement"></a>

## SlideRCGraphRingElement Objects

```python
class SlideRCGraphRingElement(RCGraphRingElement)
```

Element of `SlideRCGraphRing`.

<a id="schubmult.rings.combinatorial.slide_rc_ring.SlideRCGraphRing"></a>

## SlideRCGraphRing Objects

```python
class SlideRCGraphRing(RCGraphRing)
```

`RCGraphRing` snapped to quasi-Yamanouchi-class representatives; see the module docstring.

