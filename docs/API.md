# schubmult — package overview

This page summarizes the purpose and structure of the schubmult package and links to its main submodules and key symbols. Use the per-module pages in [docs/modules](docs/modules) and the source files for full API details.

## What schubmult is for

schubmult is a Python library for computing products, coproducts, and related operations for Schubert polynomials (ordinary, double, quantum, and variants). It integrates with sympy and provides command-line scripts for common workflows.

## Top-level package

- Source: [src/schubmult/__init__.py](src/schubmult/__init__.py)

## Mult (multiplication kernels)

Contains implementations of multiplication algorithms for different polynomial flavors.

- Package: [src/schubmult/mult](src/schubmult/mult)
- Ordinary single-variable algorithms: [`src.schubmult.mult.single.single_variable`](src/schubmult/mult/single.py), [`src.schubmult.mult.single.schubmult_py`](src/schubmult/mult/single.py)
- Double Schubert polynomials: [`src.schubmult.mult.double.mult_poly_double`](src/schubmult/mult/double.py), [`src.schubmult.mult.double.schubmult_double`](src/schubmult/mult/double.py)
- Quantum / quantum-double: [`src.schubmult.mult.quantum.schubmult_q`](src/schubmult/mult/quantum.py), [`src.schubmult.mult.quantum_double.schubmult_q_double`](src/schubmult/mult/quantum_double.py)
- Positivity helpers and transformations: [`src.schubmult.mult.positivity.posify`](src/schubmult/mult/positivity.py)

Notes:
- The `double` and `quantum_double` implementations include both core multiplication and helper routines for factorization and q-variable bookkeeping (`factor_out_q`, `mul_q_dict`).

## Rings and algebraic structures

Type-safe ring abstractions, bases, and algebra operations.

- Package: [src/schubmult/rings](src/schubmult/rings)
- Abstract polynomial types and printing: [`src.schubmult.rings.abstract_schub_poly.AbstractSchubPoly`](src/schubmult/rings/abstract_schub_poly.py)
- Concrete rings and elements: [`src.schubmult.rings.schubert_ring.DoubleSchubertRing`](src/schubmult/rings/schubert_ring.py), [`src.schubmult.rings.quantum_schubert_ring.QuantumDoubleSchubertRing`](src/schubmult/rings/quantum_schubert_ring.py)
- Polynomial bases and algebra utilities: [`src.schubmult.rings.polynomial_basis.PolynomialBasis`](src/schubmult/rings/polynomial_basis.py), [`src.schubmult.rings.polynomial_algebra.PolynomialAlgebra`](src/schubmult/rings/polynomial_algebra.py)
- Free algebra and combinatorial bases: [`src.schubmult.rings.free_algebra.FreeAlgebra`](src/schubmult/rings/free_algebra.py), [`src.schubmult.rings.free_algebra_basis.SchubertBasis`](src/schubmult/rings/free_algebra_basis.py)

## Schubert combinatorics (schub_lib)

Combinatorial objects and algorithms underlying Schubert computations.

- Package: [src/schubmult/schub_lib](src/schubmult/schub_lib)
- Permutations and codes: [`src.schubmult.schub_lib.perm_lib.Permutation`](src/schubmult/schub_lib/perm_lib.py)
- RC-graphs and crystals: [`src.schubmult.schub_lib.rc_graph.RCGraph`](src/schubmult/schub_lib/rc_graph.py)
- Schubert polynomial construction and divided difference operators: [`src.schubmult.schub_lib.schub_poly.schubpoly`](src/schubmult/schub_lib/schub_poly.py)

## Symmetric polynomial utilities

Helpers for elementary and complete symmetric polynomials used in coefficient computations.

- Package: [src/schubmult/symmetric_polynomials](src/schubmult/symmetric_polynomials)
- Elementary and complete bases: [`src.schubmult.symmetric_polynomials.elem_sym.E`](src/schubmult/symmetric_polynomials/elem_sym.py), [`src.schubmult.symmetric_polynomials.complete_sym.H`](src/schubmult/symmetric_polynomials/complete_sym.py)
- Utility functions: [`src.schubmult.symmetric_polynomials.functions.genvars`](src/schubmult/symmetric_polynomials/functions.py)

## Utilities

Common helpers used across the project.

- Package: [src/schubmult/utils](src/schubmult/utils)
- Argument parsing for scripts: [`src.schubmult.utils.argparse.schub_argparse`](src/schubmult/utils/argparse.py)
- Multiplication dict helpers: [`src.schubmult.utils._mul_utils._mul_schub_dicts`](src/schubmult/utils/_mul_utils.py)
- Permutation utilities: [`src.schubmult.utils.perm_utils.getpermval`](src/schubmult/utils/perm_utils.py)
- Logging helpers: [`src.schubmult.utils.logging.get_logger`](src/schubmult/utils/logging.py)
- Test helpers: [`src.schubmult.utils.test_utils.generate_all`](src/schubmult/utils/test_utils.py)

## Scripts / CLI

Command-line entry points for common tasks and experiments.

- Script entry points (console scripts defined in pyproject.toml):
  - `schubmult_py` → [src/scripts/schubmult_py.py](src/scripts/schubmult_py.py) (`main`)
  - `schubmult_double` → [src/scripts/schubmult_double.py](src/scripts/schubmult_double.py) (`main`)
  - `schubmult_q` → [src/scripts/schubmult_q.py](src/scripts/schubmult_q.py) (`main`)
- Script helpers and full displays: see [`src/scripts/schubmult_double._display`](src/scripts/schubmult_double.py) and [`src/scripts/schubmult_q._display_full`](src/scripts/schubmult_q.py)

## Documentation generation

Minimal generator that produces the API index and per-module files.

- Script: [`src.schubmult.generate_docs.main`](src/schubmult/generate_docs.py)
- Per-module files live in: [docs/modules](docs/modules)

## How to extend these docs

- The per-module stubs in [docs/modules](docs/modules) are generated from source and list all top-level classes and functions. To add prose descriptions, edit or replace the corresponding file in [docs/modules] with a richer markdown page that links to the source files and symbols above.
- To re-generate the index and stubs, run the generator:
  ```sh
  python3 src/schubmult/generate_docs.py
  ```
