from pathlib import Path

OUT = Path("docs/modules")
OUT.mkdir(parents=True, exist_ok=True)

files = {
    "mult_positivity.md": """# schubmult.mult.positivity — positivity and decomposition helpers

Source: src/schubmult/mult/positivity.py

Purpose
- Decompose polynomial expressions into manifestly positive combinations of basis elements and attempt algorithmic positivity proofs for coefficients arising in Schubert products.
- Used by double and quantum-double multiplication paths to produce "display-positive" output for scripts.

Primary functions and their roles
- compute_positive_rep — Set up and solve integer/linear programs (via PuLP) that split an expression into positive base vectors and a remainder; returns explicit decomposition when successful.
- posify — High-level canonicalizer transforming possibly-signed polynomial expressions into sums of positive pieces using combinatorial splitting and compute_positive_rep.
- posify_generic_partial — Adapter used by precomputed generic partial products to ensure local positivity invariants.
- dualpieri — Dual Pieri–type enumeration used internally for certain coefficient expansions.
- Misc helpers: shiftsub, shiftsubz, split_flat_term, flatten_factors, find_base_vectors, etc.

Notes for maintainers
- Routines mix symbolic rearrangement and integer-program solving; see compute_positive_rep for solver lifecycle and solver choices.
- posify is cached and tuned for common Schubert sizes; altering canonicalization or variable ordering can affect cache behavior.
""",
    "mult_double.md": """# schubmult.mult.double — double Schubert multiplication kernels

Source: src/schubmult/mult/double.py

Purpose
- Algorithms to compute products and coproducts for double Schubert polynomials and ancillary helpers used by higher-level ring classes.

Key routines
- schubmult_double — Core multiplication entry point for double Schubert polynomial products (permutation × permutation); returns coefficient dictionaries keyed by permutations.
- schubmult_double_pair — Pairwise multiplication routine used for basis elements.
- schub_coprod_double — Compute permutation-level coproduct for double Schubert polynomials.
- Helpers: mult_poly_double, single_variable, nilhecke_mult, and generic enumerators (schubmult_double_pair_generic, _alt variants).

Integration notes
- These kernels call into positivity helpers for display-positive output.
- Functions accept generating-set / variable arguments to control double-variable output and substitution behavior.
""",
    "perm_lib.md": """# schubmult.schub_lib.perm_lib — permutation primitives

Source: src/schubmult/schub_lib/perm_lib.py

Overview
- Lightweight Permutation class and utilities for Lehmer/code conversions, Bruhat relations, cycles, and coset representatives used across the codebase.

Primary symbol
- Permutation — core permutation type with helpers:
  - construction: from code/cycles, uncode, code
  - properties: inv, code/trimcode, descents, cycles
  - comparators and Bruhat utilities used by algorithms that enumerate or compare permutations.

Useful utilities
- code, uncode, permtrim, theta, strict_theta, longest_element — common transformations.
- cyclic_sort, cyclic_sort_min, get_cycles — cycle and ordering helpers.

Notes
- Permutation objects are treated as immutable and often cached for allocation savings; many algorithms rely on fast equality and code-based comparisons.
""",
    "rc_graph.md": """# schubmult.schub_lib.rc_graph — RC-graphs and crystal interactions

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
""",
    "schubert_ring.md": """# schubmult.rings.schubert_ring — Schubert ring classes and elements

Source: src/schubmult/rings/schubert_ring.py

Purpose
- Concrete ring classes exposing Schubert-basis algebras (single and double flavors) with dict-backed elements supporting arithmetic, pretty printing, and integrations with multiplication kernels.

Primary types
- DoubleSchubertRing — ring container that constructs DoubleSchubertElement instances and wires multiplication to kernels in mult.double.
- DoubleSchubertElement — dict-backed element representing Schubert-basis expressions; supports addition, multiplication, substitution, and canonical rendering.

Helpers and factories
- DSx — factory for basis elements.
- Integration points to elem_sym and polynomial basis utilities for coefficient expansions.

Notes
- Elements inherit from BaseSchubertElement; equality/approximate comparisons account for symbolic coefficient expansions and substitution rules.
""",
    "perm_utils.md": """# schubmult.utils.perm_utils — permutation helpers and combinatorial helpers

Source: src/schubmult/utils/perm_utils.py

Purpose
- Helper functions for permutation-based combinatorics: Bruhat tests, artin sequences, parabolic detection, dictionary merges and small counting utilities used throughout mult and schub_lib.

Key functions
- getpermval — safe retrieval from permutation-like sequences.
- add_perm_dict — merge/add two permutation-keyed coefficient dictionaries.
- count_less_than — utility used in inversion and counting computations.
- has_bruhat_descent / has_bruhat_ascent / count_bruhat / is_parabolic — combinatorial predicates.
- omega, sg, p_trans, mu_A, get_cycles — transformation and cycle utilities.

Notes
- Many algorithms rely on these predicates for branching and enumerations; keep edge-case tests for small/empty permutations and short sequences.
""",
    "schubmult_double_script.md": """# scripts.schubmult_double — CLI for double Schubert products

Source: src/scripts/schubmult_double.py

Purpose
- Command-line tool for computing and pretty-printing double Schubert polynomial products; supports display-positive mode via posify.

Important functions
- _display_full — format and print full coefficient dictionaries; prepare substitution maps for nicer variable rendering.
- pre_posify / sv_posify — wrappers that convert raw symbolic coefficients into display-positive canonical forms by delegating to mult.positivity.posify.
- main — CLI entrypoint wired by pyproject scripts (schubmult_double).

Usage notes
- The script supports flags for coproducts, display-positive output, and variable control. See tests in src/tests/script_tests for example invocations.
""",
}

for name, content in files.items():
    p = OUT / name
    p.write_text(content, encoding="utf8")
    print("wrote", p)