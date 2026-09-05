# Changelog

## 4.2.0

### Fixed

- **Parabolic quantum ring constructors were unusable due to a circular import.**
  `schubmult.rings.schubert.schubert_ring` imported `quantum_schubert_ring` at module
  load time, closing a cycle through `parabolic_quantum_double_schubert_ring`. Entering
  the cycle at the parabolic module — as any of `ParabolicQuantumDoubleSchubertRing`,
  `make_parabolic_quantum_basis`, or `make_single_parabolic_quantum_basis` did — raised

  ```
  ImportError: cannot import name 'ParabolicQuantumDoubleSchubertElement' from
  partially initialized module '...parabolic_quantum_double_schubert_ring'
  ```

  The import is now deferred to the single call site that needs it. Reported by an
  external smoke test of 4.1.0, which found the bug while trying to compute QH*(P^3);
  that computation is now a regression test.

- **Five `dual_basis()` declarations raised instead of returning the dual basis.**
  Each imported the right class name from the wrong module:
  `GrothendieckPolyBasis` looked for `GrothendieckBasis` in `free_algebra.schubert_basis`,
  `GlidePolyBasis` and `LascouxPolyBasis` both looked in `free_algebra.fundamental_slide_basis`,
  and `LascouxBasis` looked for `LascouxPolyBasis` in `polynomial_algebra.key_poly_basis`.
  `GrovePolyBasis` declared no dual at all and inherited the base implementation, which
  returns `None`. Every declared dual pair now round-trips in both directions.

### Added

- Grothendieck polynomial multiplication (`schubmult.mult.groth`,
  `schubmult.mult.groth_double`) and a double Grothendieck ring
  (`rings.schubert.double_grothendieck_ring`).
- Lascoux and glide polynomial bases, in both the free algebra
  (`rings.free_algebra.lascoux_basis`, `glide_basis`) and the polynomial algebra
  (`rings.polynomial_algebra.lascoux_poly_basis`, `glide_poly_basis`).
- Chevalley/Monk formula module `rings.schubert.chevalley`.
- Double Grothendieck polynomial separated-descents multiplication (`schubmult.mult.separated_descents`).
- Weak order operations on `Permutation`: `weak_order_leq`, `weak_order_meet`,
  `weak_order_join`.
- Verification scripts `groth_lr_rule` and `verify_quantum_triple_positive`.

### Removed

- `rings.combinatorial.bounded_double_rc_factor_algebra`, an abandoned line of work
  that could not be imported: it required `DecoratedRCGraph`, which had already been
  commented out of `combinatorics.rc_graph`. Its three classes
  (`BoundedDoubleRCFactorAlgebra`, `BoundedDoubleRCFactorAlgebraElement`,
  `BoundedDoubleRCFactorPrintingTerm`) were advertised by `dir(schubmult)` but raised
  `AttributeError` on access. Nothing referenced them. The unrelated
  `BoundedRCFactorAlgebra` is unaffected.
- The commented-out `DecoratedRCGraph` class body and its `full_CEM_double` /
  `elem_squash` helpers, which had no live definitions or callers.

### Notes

- The free algebra is the dual of the polynomial algebra, so each free algebra basis
  is the dual basis to the correspondingly named polynomial basis. The naming
  convention carries the pairing: `SchubertBasis` / `SchubertPolyBasis`,
  `GrothendieckBasis` / `GrothendieckPolyBasis`, `LascouxBasis` / `LascouxPolyBasis`,
  `GlideBasis` / `GlidePolyBasis`, `GroveBasis` / `GrovePolyBasis`, `KeyBasis` /
  `KeyPolyBasis`, `ForestBasis` / `ForestPolyBasis`, `FundamentalSlideBasis` /
  `FundamentalSlidePolyBasis`, `MonomialSlideBasis` / `MonomialSlidePolyBasis`,
  `CompositionSchubertBasis` / `CompositionSchubertPolyBasis`,
  `SeparatedDescentsBasis` / `SepDescPolyBasis`, `ElementaryBasis` /
  `ElemSymPolyBasis`, and `WordBasis` / `MonomialBasis`. Where a dual is itself
  realized inside the free algebra it is named explicitly: `ForestDual`, `GlideDual`,
  `GroveDual`.
- The PyPI project description is generated from `README.md`, which documents the
  current (4.x) API. Releases before this one carried a v1.x-era description.
