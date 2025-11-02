# schubmult.mult.positivity — positivity and decomposition helpers

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
