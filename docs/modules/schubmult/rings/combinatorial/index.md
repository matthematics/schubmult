<a id="schubmult.rings.combinatorial"></a>

# schubmult.rings.combinatorial

Combinatorial rings: rings whose basis elements are combinatorial objects (RC graphs, BPDs, WC graphs,
tableaux, ...) rather than permutations.

The central object is `RCGraphRing`; most other rings here are quotients or variants of it that snap
products to canonical representatives (`HWRCGraphRing`, `KeyRCGraphRing`, `QYRCGraphRing`,
`ForestRCGraphRing`, `SlideRCGraphRing`, ...). `BoundedRCFactorAlgebra` and `GrassTensorAlgebra`
provide factorizations into Grassmannian pieces used to compute products and coproducts.

