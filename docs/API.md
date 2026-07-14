# schubmult module guide

This page is the friendly map of the package.

If you want a place to start as a user, read this page. If you want the fuller callable-level reference, use [PACKAGE_API.md](PACKAGE_API.md). If you want generated per-module stubs, browse [docs/modules](modules).

## The short version

Most users only need four parts of the package:

- [src/schubmult/__init__.py](../src/schubmult/__init__.py) for top-level imports like `Permutation`, `Sx`, `DSx`, `RCGraph`, and `RootTableau`
- [src/schubmult/rings/schubert](../src/schubmult/rings/schubert) for ordinary, double, Grothendieck, and quantum Schubert rings
- [src/schubmult/combinatorics](../src/schubmult/combinatorics) for permutation and pipe-dream style objects
- [src/schubmult/_scripts](../src/schubmult/_scripts) for command-line tools and verification scripts

If you are unsure where to start, start with the ring interface, not the low-level multiplication kernels.

## Which module should I use?

| If you want to... | Start here | Why |
|---|---|---|
| Multiply Schubert classes as a user | [src/schubmult/rings/schubert](../src/schubmult/rings/schubert) | This is the highest-level algebra interface |
| Work directly with permutations | [src/schubmult/combinatorics/permutation.py](../src/schubmult/combinatorics/permutation.py) | `Permutation` is the basic indexing object used everywhere |
| Explore RC-graphs, pipe dreams, tableaux, or crystals | [src/schubmult/combinatorics](../src/schubmult/combinatorics) | This is the main combinatorics layer |
| Call the fast multiplication algorithms directly | [src/schubmult/mult](../src/schubmult/mult) | These are the kernel implementations under the ring wrappers |
| Build or compare polynomial bases | [src/schubmult/rings/polynomial_algebra](../src/schubmult/rings/polynomial_algebra) | Basis-oriented polynomial algebra lives here |
| Work with noncommutative bases and free algebra constructions | [src/schubmult/rings/free_algebra](../src/schubmult/rings/free_algebra) | This is the free-algebra side of the project |
| Use combinatorial objects as algebra elements | [src/schubmult/rings/combinatorial](../src/schubmult/rings/combinatorial) | Rings built from RC-graphs, BPDs, plactic data, and related models |
| Use the package from a shell | [src/schubmult/_scripts](../src/schubmult/_scripts) | Console entry points and research scripts live here |
| Stay compatible with SymPy and SymEngine | [src/schubmult/symbolic](../src/schubmult/symbolic) | Shared symbolic helpers and variable factories |

## Recommended entry points

### 1. Top-level imports

[src/schubmult/__init__.py](../src/schubmult/__init__.py) lazily re-exports the objects most users care about.

Use this when you want a convenient import style such as:

```python
from schubmult import Permutation, Sx, DSx, RCGraph, RootTableau
```

This is the right starting point for interactive use, notebooks, and quick experiments.

### 2. Schubert ring modules

[src/schubmult/rings/schubert](../src/schubmult/rings/schubert) is the main user-facing algebra layer.

Important modules:

- [src/schubmult/rings/schubert/schubert_ring.py](../src/schubmult/rings/schubert/schubert_ring.py): ordinary and double Schubert rings used through `Sx` and `DSx`
- [src/schubmult/rings/schubert/double_schubert_ring.py](../src/schubmult/rings/schubert/double_schubert_ring.py): dedicated double Schubert ring implementation details
- [src/schubmult/rings/schubert/grothendieck_ring.py](../src/schubmult/rings/schubert/grothendieck_ring.py): Grothendieck ring functionality
- [src/schubmult/rings/schubert/quantum_schubert_ring.py](../src/schubmult/rings/schubert/quantum_schubert_ring.py): quantum Schubert rings
- [src/schubmult/rings/schubert/quantum_double_schubert_ring.py](../src/schubmult/rings/schubert/quantum_double_schubert_ring.py): quantum double rings
- [src/schubmult/rings/schubert/parabolic_quantum_schubert_ring.py](../src/schubmult/rings/schubert/parabolic_quantum_schubert_ring.py): parabolic quantum variants
- [src/schubmult/rings/schubert/parabolic_quantum_double_schubert_ring.py](../src/schubmult/rings/schubert/parabolic_quantum_double_schubert_ring.py): parabolic quantum double variants
- [src/schubmult/rings/schubert/base_schubert_ring.py](../src/schubmult/rings/schubert/base_schubert_ring.py): shared infrastructure used by the Schubert-family rings

Use this part of the package when your mental model is “I want ring elements that I can add, multiply, and expand.”

### 3. Combinatorics modules

[src/schubmult/combinatorics](../src/schubmult/combinatorics) is where the underlying discrete objects live.

Best starting points:

- [src/schubmult/combinatorics/permutation.py](../src/schubmult/combinatorics/permutation.py): the core indexing object for Schubert data
- [src/schubmult/combinatorics/rc_graph.py](../src/schubmult/combinatorics/rc_graph.py): reduced compatible graphs and crystal-style operations
- [src/schubmult/combinatorics/root_tableau.py](../src/schubmult/combinatorics/root_tableau.py): root tableaux and dual Knuth style structures
- [src/schubmult/combinatorics/bpd.py](../src/schubmult/combinatorics/bpd.py): bumpless pipe dreams
- [src/schubmult/combinatorics/hpd.py](../src/schubmult/combinatorics/hpd.py): hybrid pipe dreams
- [src/schubmult/combinatorics/pipe_dream.py](../src/schubmult/combinatorics/pipe_dream.py): classical pipe dream representation
- [src/schubmult/combinatorics/crystal_graph.py](../src/schubmult/combinatorics/crystal_graph.py): crystal interface and tensor-style behavior
- [src/schubmult/combinatorics/wc_graph.py](../src/schubmult/combinatorics/wc_graph.py): weak-compatible graph structures

Other modules in this package are useful, but many are more specialized: insertion models, tableau variants, alternative graph models, and research-specific combinatorial objects.

### 4. Multiplication kernels

[src/schubmult/mult](../src/schubmult/mult) contains the lower-level algorithms used by the ring classes.

Important modules:

- [src/schubmult/mult/single.py](../src/schubmult/mult/single.py): ordinary Schubert multiplication
- [src/schubmult/mult/double.py](../src/schubmult/mult/double.py): double and mixed-variable multiplication
- [src/schubmult/mult/quantum.py](../src/schubmult/mult/quantum.py): quantum multiplication
- [src/schubmult/mult/quantum_double.py](../src/schubmult/mult/quantum_double.py): quantum double multiplication
- [src/schubmult/mult/positivity.py](../src/schubmult/mult/positivity.py): positivity display helpers and related utilities

Use these modules when you want direct algorithm access, performance-focused experiments, or implementation details behind ring multiplication.

## The ring subpackages

The [src/schubmult/rings](../src/schubmult/rings) package is broader than just Schubert rings.

### `rings.schubert`

The main place for Schubert-family rings. This is the best entry point for most users.

### `rings.polynomial_algebra`

[src/schubmult/rings/polynomial_algebra](../src/schubmult/rings/polynomial_algebra) is for basis-driven polynomial algebra.

Good modules to know:

- [src/schubmult/rings/polynomial_algebra/base_polynomial_basis.py](../src/schubmult/rings/polynomial_algebra/base_polynomial_basis.py): shared basis protocol
- [src/schubmult/rings/polynomial_algebra/schubert_poly_basis.py](../src/schubmult/rings/polynomial_algebra/schubert_poly_basis.py): Schubert basis
- [src/schubmult/rings/polynomial_algebra/monomial_basis.py](../src/schubmult/rings/polynomial_algebra/monomial_basis.py): monomial basis
- [src/schubmult/rings/polynomial_algebra/elem_sym_poly_basis.py](../src/schubmult/rings/polynomial_algebra/elem_sym_poly_basis.py): elementary symmetric basis
- [src/schubmult/rings/polynomial_algebra/key_poly_basis.py](../src/schubmult/rings/polynomial_algebra/key_poly_basis.py): key polynomials
- [src/schubmult/rings/polynomial_algebra/sepdesc_poly_basis.py](../src/schubmult/rings/polynomial_algebra/sepdesc_poly_basis.py): separated-descents basis

This is the part of the package to use when your question is about changing bases, expanding in a basis, or building new polynomial-basis objects.

### `rings.free_algebra`

[src/schubmult/rings/free_algebra](../src/schubmult/rings/free_algebra) contains noncommutative basis families and free-algebra constructions.

Representative modules:

- [src/schubmult/rings/free_algebra/free_algebra_basis.py](../src/schubmult/rings/free_algebra/free_algebra_basis.py)
- [src/schubmult/rings/free_algebra/word_basis.py](../src/schubmult/rings/free_algebra/word_basis.py)
- [src/schubmult/rings/free_algebra/schubert_basis.py](../src/schubmult/rings/free_algebra/schubert_basis.py)
- [src/schubmult/rings/free_algebra/grothendieck_basis.py](../src/schubmult/rings/free_algebra/grothendieck_basis.py)
- [src/schubmult/rings/free_algebra/j_basis.py](../src/schubmult/rings/free_algebra/j_basis.py)
- [src/schubmult/rings/free_algebra/jt_basis.py](../src/schubmult/rings/free_algebra/jt_basis.py)

Use this when you want free-algebra bases rather than commutative polynomial or Schubert ring elements.

### `rings.combinatorial`

[src/schubmult/rings/combinatorial](../src/schubmult/rings/combinatorial) turns combinatorial objects into algebra elements.

Common examples:

- [src/schubmult/rings/combinatorial/rc_graph_ring.py](../src/schubmult/rings/combinatorial/rc_graph_ring.py)
- [src/schubmult/rings/combinatorial/wc_graph_ring.py](../src/schubmult/rings/combinatorial/wc_graph_ring.py)
- [src/schubmult/rings/combinatorial/bpd_ring.py](../src/schubmult/rings/combinatorial/bpd_ring.py)
- [src/schubmult/rings/combinatorial/plactic_algebra.py](../src/schubmult/rings/combinatorial/plactic_algebra.py)
- [src/schubmult/rings/combinatorial/dual_rc_graph_ring.py](../src/schubmult/rings/combinatorial/dual_rc_graph_ring.py)
- [src/schubmult/rings/combinatorial/bounded_rc_factor_algebra.py](../src/schubmult/rings/combinatorial/bounded_rc_factor_algebra.py)

This subpackage is especially useful for research workflows where the combinatorial model itself is the object of study.

### Other ring-level modules

The remaining modules in [src/schubmult/rings](../src/schubmult/rings) provide general algebra infrastructure:

- [src/schubmult/rings/base_ring.py](../src/schubmult/rings/base_ring.py)
- [src/schubmult/rings/tensor_ring.py](../src/schubmult/rings/tensor_ring.py)
- [src/schubmult/rings/direct_product_ring.py](../src/schubmult/rings/direct_product_ring.py)
- [src/schubmult/rings/product_ring.py](../src/schubmult/rings/product_ring.py)
- [src/schubmult/rings/nsym.py](../src/schubmult/rings/nsym.py)
- [src/schubmult/rings/quasisymmetric_functions.py](../src/schubmult/rings/quasisymmetric_functions.py)
- [src/schubmult/rings/thompson_algebra.py](../src/schubmult/rings/thompson_algebra.py)

Most users will only visit these after they already know which higher-level ring family they want.

## Symbolic and utility modules

### `symbolic`

[src/schubmult/symbolic](../src/schubmult/symbolic) is the compatibility layer for symbolic arithmetic.

Use this package when you need variables, symbolic expressions, or helpers that should work consistently with the rest of `schubmult`.

Helpful modules:

- [src/schubmult/symbolic/__init__.py](../src/schubmult/symbolic/__init__.py): shared imports and symbolic facade
- [src/schubmult/symbolic/common_polys.py](../src/schubmult/symbolic/common_polys.py): common polynomial helpers
- [src/schubmult/symbolic/poly/variables.py](../src/schubmult/symbolic/poly/variables.py): variable factories such as `GeneratingSet`
- [src/schubmult/symbolic/symmetric_polynomials](../src/schubmult/symbolic/symmetric_polynomials): elementary, complete, and related symmetric-polynomial code

### `utils`

[src/schubmult/utils](../src/schubmult/utils) contains support code that other modules lean on.

Useful modules include:

- [src/schubmult/utils/argparse.py](../src/schubmult/utils/argparse.py): shared CLI parsing
- [src/schubmult/utils/parsing.py](../src/schubmult/utils/parsing.py): parsing helpers
- [src/schubmult/utils/perm_utils.py](../src/schubmult/utils/perm_utils.py): permutation and coefficient-dictionary utilities
- [src/schubmult/utils/schub_lib.py](../src/schubmult/utils/schub_lib.py): legacy helper functions used across the package
- [src/schubmult/utils/tuple_utils.py](../src/schubmult/utils/tuple_utils.py): tuple-level helpers

Most end users should treat `utils` as internal support code unless a specific helper has already been recommended in examples or docs.

## Scripts and command-line tools

[src/schubmult/_scripts](../src/schubmult/_scripts) serves two different roles:

- the actual console entry points installed with the package
- a large collection of exploratory and verification scripts used during development and research

Installed entry-point style scripts include:

- [src/schubmult/_scripts/schubmult_py.py](../src/schubmult/_scripts/schubmult_py.py)
- [src/schubmult/_scripts/grothmult_py.py](../src/schubmult/_scripts/grothmult_py.py)
- [src/schubmult/_scripts/schubmult_double.py](../src/schubmult/_scripts/schubmult_double.py)
- [src/schubmult/_scripts/schubmult_q.py](../src/schubmult/_scripts/schubmult_q.py)
- [src/schubmult/_scripts/schubmult_q_double.py](../src/schubmult/_scripts/schubmult_q_double.py)
- [src/schubmult/_scripts/schubert_formula.py](../src/schubmult/_scripts/schubert_formula.py)
- [src/schubmult/_scripts/schubprompt.py](../src/schubmult/_scripts/schubprompt.py)
- [src/schubmult/_scripts/lr_rule_verify.py](../src/schubmult/_scripts/lr_rule_verify.py)

Many other files in this directory are research scripts rather than polished public APIs. They can still be valuable examples, but they are not the best first stop for learning the package.

## Stable starting points versus deeper internals

If you are documenting or teaching the package, treat these as the stable starting points:

- [src/schubmult/__init__.py](../src/schubmult/__init__.py)
- [src/schubmult/combinatorics/permutation.py](../src/schubmult/combinatorics/permutation.py)
- [src/schubmult/rings/schubert](../src/schubmult/rings/schubert)
- [src/schubmult/mult](../src/schubmult/mult)
- [src/schubmult/symbolic](../src/schubmult/symbolic)

Treat these as more specialized or research-oriented layers:

- [src/schubmult/rings/combinatorial](../src/schubmult/rings/combinatorial)
- [src/schubmult/rings/free_algebra](../src/schubmult/rings/free_algebra)
- most of [src/schubmult/_scripts](../src/schubmult/_scripts)
- specialized objects in [src/schubmult/combinatorics](../src/schubmult/combinatorics) beyond permutations, RC-graphs, root tableaux, and pipe dreams

## Suggested learning order

If you want to learn the codebase without getting lost, this order works well:

1. [src/schubmult/__init__.py](../src/schubmult/__init__.py)
2. [src/schubmult/combinatorics/permutation.py](../src/schubmult/combinatorics/permutation.py)
3. [src/schubmult/rings/schubert/schubert_ring.py](../src/schubmult/rings/schubert/schubert_ring.py)
4. [src/schubmult/mult/single.py](../src/schubmult/mult/single.py)
5. [src/schubmult/combinatorics/rc_graph.py](../src/schubmult/combinatorics/rc_graph.py)
6. [src/schubmult/combinatorics/root_tableau.py](../src/schubmult/combinatorics/root_tableau.py)
7. whichever specialized ring or combinatorics subpackage matches your problem

## Exhaustive references

For more detail after this overview:

- [PACKAGE_API.md](PACKAGE_API.md) for a fuller practical reference
- [docs/modules](modules) for generated per-module pages

To regenerate the module stubs:

```sh
python3 src/schubmult/generate_docs.py
```
