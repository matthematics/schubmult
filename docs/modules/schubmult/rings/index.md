<a id="schubmult.rings"></a>

# schubmult.rings

Ring structures built on combinatorial bases.

Subpackages:

- `schubert`: the Schubert-family rings (``Sx``, ``DSx``, ``Gx``, ``QSx``, ...) -- the
  main user-facing algebra interface.
- `polynomial_algebra`: the polynomial ring ``Z[x_1, x_2, ...]`` with pluggable bases
  (monomial, Schubert, key, slide, Grothendieck, ...).
- `free_algebra`: noncommutative free algebras with combinatorial bases.
- `combinatorial`: rings whose basis elements are combinatorial objects (RC graphs,
  BPDs, plactic classes, ...) rather than permutations.

Modules here: `base_ring` (the shared dict-based ring/element machinery), `tensor_ring`,
`product_ring`, `direct_product_ring`, `printing` (display symbols), `nsym`,
`quasisymmetric_functions`, `thompson_algebra`.

This ``__init__`` re-exports nothing (the commented-out block below is legacy);
import from the subpackages directly.

