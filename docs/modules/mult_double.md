# schubmult.mult.double — double Schubert multiplication kernels

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
