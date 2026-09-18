<a id="schubmult.rings.combinatorial.schubert_rc_ring"></a>

# schubmult.rings.combinatorial.schubert\_rc\_ring

`SchubertRCGraphRing`: `RCGraphRing` whose product is computed through `BoundedRCFactorAlgebra`.

Each RC graph is factored into elementary-symmetric RC graphs (``_factor_rc``), the factors
are multiplied in the bounded factor algebra, and the result is converted back; ``schubert_poly``
is the sum of all RC graphs of a permutation.

<a id="schubmult.rings.combinatorial.schubert_rc_ring.SchubertRCGraphRingElement"></a>

## SchubertRCGraphRingElement Objects

```python
class SchubertRCGraphRingElement(RCGraphRingElement)
```

Element of `SchubertRCGraphRing`.

<a id="schubmult.rings.combinatorial.schubert_rc_ring.SchubertRCGraphRing"></a>

## SchubertRCGraphRing Objects

```python
class SchubertRCGraphRing(RCGraphRing)
```

`RCGraphRing` multiplying via `BoundedRCFactorAlgebra` factorizations; see the module docstring.

