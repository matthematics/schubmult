# schubmult.utils.perm_utils — permutation helpers and combinatorial helpers

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
