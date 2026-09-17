<a id="schubmult.rings.polynomial_algebra.composition_schubert_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.composition\_schubert\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.composition_schubert_poly_basis.CompositionSchubertPolyBasis"></a>

## CompositionSchubertPolyBasis Objects

```python
class CompositionSchubertPolyBasis(SchubertPolyBasis)
```

Wrapper basis for Schubert polynomials indexed by weak compositions.

Keys are weak compositions interpreted as Lehmer codes. The underlying
computations are delegated to :class:`SchubertPolyBasis`, while display
uses the same Schubert printing as the corresponding permutation key.

