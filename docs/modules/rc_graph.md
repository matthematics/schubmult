# schubmult.schub_lib.rc_graph — RC-graphs and crystal interactions

Source: src/schubmult/schub_lib/rc_graph.py

Purpose
- RCGraph models reduced-compatible graphs underlying Schubert combinatorics; provides crystal operators, tableau decompositions and combinatorial multiplication (Monk-type) in the crystal language.

Key class and methods
- RCGraph
  - Crystal API: raising/lowering operators, epsilon/phi, crystal_weight, to_lowest_weight, is_highest_weight.
  - Combinatorial products: monk_crystal_mul for Monk product.
  - Structural helpers: tableau_decomp, all_rc_graphs, rowrange, vertical_cut.
  - Printing: human-readable grid representations for debugging and tests.

Notes
- RCGraph integrates with FreeAlgebra and NilHecke components; check the test suite and scripts (e.g., scripts/monk_crystal.py) for typical usage patterns.
