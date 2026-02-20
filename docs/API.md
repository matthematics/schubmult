# schubmult — API and learning map

This page gives a practical map of where things live in the current codebase, plus the best order to learn them.

## Start here (new learners)

If you are learning Schubert calculus from this repo, use this order:

1. [Beginner guide](learn_schubert_calculus.md)
2. [Practice exercises](learn_schubert_calculus_exercises.md)
3. This API map (to connect ideas to source code)

## Top-level package

- Main package: [src/schubmult](../src/schubmult)
- Package init: [src/schubmult/__init__.py](../src/schubmult/__init__.py)

## Combinatorics layer

Core combinatorial objects used throughout the project.

- Package: [src/schubmult/combinatorics](../src/schubmult/combinatorics)
- Permutations: [src/schubmult/combinatorics/permutation.py](../src/schubmult/combinatorics/permutation.py)
- RC-graphs: [src/schubmult/combinatorics/rc_graph.py](../src/schubmult/combinatorics/rc_graph.py)
- Root tableaux: [src/schubmult/combinatorics/root_tableau.py](../src/schubmult/combinatorics/root_tableau.py)
- Schubert polynomial helpers: [src/schubmult/combinatorics/schub_poly.py](../src/schubmult/combinatorics/schub_poly.py)

## Multiplication kernels

Algorithmic multiplication routines for different variants.

- Package: [src/schubmult/mult](../src/schubmult/mult)
- Ordinary: [src/schubmult/mult/single.py](../src/schubmult/mult/single.py)
- Double: [src/schubmult/mult/double.py](../src/schubmult/mult/double.py)
- Quantum: [src/schubmult/mult/quantum.py](../src/schubmult/mult/quantum.py)
- Quantum double: [src/schubmult/mult/quantum_double.py](../src/schubmult/mult/quantum_double.py)
- Positivity display utilities: [src/schubmult/mult/positivity.py](../src/schubmult/mult/positivity.py)

## Ring implementations

Concrete algebraic ring types and their elements.

- Package: [src/schubmult/rings/schubert](../src/schubmult/rings/schubert)
- Single/ordinary ring: [src/schubmult/rings/schubert/schubert_ring.py](../src/schubmult/rings/schubert/schubert_ring.py)
- Double ring: [src/schubmult/rings/schubert/double_schubert_ring.py](../src/schubmult/rings/schubert/double_schubert_ring.py)
- Quantum ring (compat entry): [src/schubmult/rings/schubert/quantum_schubert_ring.py](../src/schubmult/rings/schubert/quantum_schubert_ring.py)
- Quantum double ring: [src/schubmult/rings/schubert/quantum_double_schubert_ring.py](../src/schubmult/rings/schubert/quantum_double_schubert_ring.py)
- Parabolic quantum ring: [src/schubmult/rings/schubert/parabolic_quantum_schubert_ring.py](../src/schubmult/rings/schubert/parabolic_quantum_schubert_ring.py)
- Parabolic quantum double ring: [src/schubmult/rings/schubert/parabolic_quantum_double_schubert_ring.py](../src/schubmult/rings/schubert/parabolic_quantum_double_schubert_ring.py)

## Polynomial algebra and bases

Abstract and concrete polynomial basis classes.

- Package: [src/schubmult/rings/polynomial_algebra](../src/schubmult/rings/polynomial_algebra)
- Abstract base class: [src/schubmult/rings/polynomial_algebra/base_polynomial_basis.py](../src/schubmult/rings/polynomial_algebra/base_polynomial_basis.py)
- Monomial basis: [src/schubmult/rings/polynomial_algebra/monomial_basis.py](../src/schubmult/rings/polynomial_algebra/monomial_basis.py)
- Schubert basis: [src/schubmult/rings/polynomial_algebra/schubert_poly_basis.py](../src/schubmult/rings/polynomial_algebra/schubert_poly_basis.py)
- Separated-descents basis: [src/schubmult/rings/polynomial_algebra/sepdesc_poly_basis.py](../src/schubmult/rings/polynomial_algebra/sepdesc_poly_basis.py)
- Elementary-symmetric basis: [src/schubmult/rings/polynomial_algebra/elem_sym_poly_basis.py](../src/schubmult/rings/polynomial_algebra/elem_sym_poly_basis.py)
- Compatibility exports: [src/schubmult/rings/polynomial_algebra/polynomial_basis.py](../src/schubmult/rings/polynomial_algebra/polynomial_basis.py)

## Symbolic layer

Symbolic helpers and symmetric polynomial utilities.

- Package: [src/schubmult/symbolic](../src/schubmult/symbolic)
- Schubert polynomial symbolic helpers: [src/schubmult/symbolic/poly/schub_poly.py](../src/schubmult/symbolic/poly/schub_poly.py)
- Symmetric polynomial package: [src/schubmult/symbolic/symmetric_polynomials](../src/schubmult/symbolic/symmetric_polynomials)

## CLI scripts

Console tools for computation and verification.

- Script package: [src/schubmult/_scripts](../src/schubmult/_scripts)
- Ordinary multiplication script: [src/schubmult/_scripts/schubmult_py.py](../src/schubmult/_scripts/schubmult_py.py)
- Double multiplication script: [src/schubmult/_scripts/schubmult_double.py](../src/schubmult/_scripts/schubmult_double.py)
- Quantum multiplication script: [src/schubmult/_scripts/schubmult_q.py](../src/schubmult/_scripts/schubmult_q.py)
- Quantum double script: [src/schubmult/_scripts/schubmult_q_double.py](../src/schubmult/_scripts/schubmult_q_double.py)
- Verification script: [src/schubmult/_scripts/lr_rule_verify.py](../src/schubmult/_scripts/lr_rule_verify.py)

## Regenerating API stubs

Generated per-module docs live in [docs/modules](modules).

To regenerate docs stubs:

```sh
python3 src/schubmult/generate_docs.py
```
