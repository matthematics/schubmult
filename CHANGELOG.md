# Changelog

## Unreleased

### Changed

- The PuLP requirement is now `PuLP[cbc]>=3.3.2`. `--display-positive` builds variables
  with `LpProblem.add_variable` and solves with `COIN_CMD` (the CBC binary comes from the
  `[cbc]` extra), replacing `LpVariable(...)` and `PULP_CBC_CMD`, which PuLP deprecates
  ahead of 4.0.

## 5.1.0

Grothendieck products on the command line and in `schubmult.mult`, including the quantum
and quantum double cases; the multiplication kernels stop expanding coefficients.

### Added

- **`grothmult_py`, `grothmult_double`, `grothmult_q`, `grothmult_q_double` command-line
  scripts.** The README and API reference on GitHub have described the Grothendieck CLI
  since 5.0.0 was tagged, but the 5.0.0 release only shipped `grothmult_py` (via the slow
  ring-expansion route) and `grothmult_double`; `grothmult_q` and `grothmult_q_double` did
  not exist. All four are now installed alongside the `schubmult_*` scripts, accept the
  same options (`--code`, `--display-mode`, `--mixed-var`, `--display-positive` for the
  double case), and have script tests with JSON fixtures.
- **`schubmult.mult.groth.grothmult_py(perm_dict, v, beta)`**: ordinary
  $\beta$-Grothendieck products in the $G$ basis, computed by expanding $G_v$ into Schubert
  polynomials and pushing each through the theta-code v-path layers with the binomial
  closed form `groth_elem_sym_coeff`. Replaces the ring-level `groth_mul_full_with_ring`
  route in `GrothendieckRing`; typical products are orders of magnitude faster.
- **`schubmult.mult.groth_quantum.grothmult_q`** and
  **`schubmult.mult.groth_quantum_double.grothmult_q_double`**: quantum and quantum double
  $\beta$-Grothendieck products for the Lenart-Maeno quantization $G^q_v = Q(G_v)$, with
  $Q_j = \beta^2 q_j$ so that $\beta = 0$ recovers the quantum double Schubert polynomials
  and $\beta = -1$ the Lenart-Maeno polynomials representing $QK_T(Fl_n)$. The rule is an
  equivariant quantum $K$-Pieri formula over Naito-Sagaki chains in the quantum Bruhat
  graph, assembled Molev-Sagan style. It is **conjectural**: verified symbolically in
  $QK_T(Fl_3)$ and under random specializations in $QK_T(Fl_4)$ against the
  Maeno-Naito-Sagaki presentation, and its specializations ($q = 0$, $\beta = 0$,
  $y = 0$ with $\beta = -1$) are theorems or existing kernels. Also exported:
  `grothmult_q_pieri`, `grothmult_q_double_pieri`, `grothmult_q_double_top`,
  `lm_quantize`, `qgroth_poly`, `quantum_elem_sym`, `quantum_pieri_chains`.
- `DoubleGrothendieckElement.simplify()` puts the rational coefficients in $y$ and $\beta$
  in normal form (`sympy.cancel`, then `factor`), dropping terms that simplify to zero.
  `BaseRingElement.simplify()` applies `.simplify()` coefficientwise for any ring.
- `NilHeckeRing.isobaric(..., neg=True)` and `g_isobaric(neg=True)` for the
  $\partial_i(1 - \beta x_{i+1})$ convention.
- Test coverage for every kernel in `schubmult.mult` (`tests/mult/`), the C++
  acceleration layer (`test_accel.py`), and positivity (`test_positivity.py`).
- API reference regenerated from docstrings with pydoc-markdown and published with
  MkDocs (`docs/`, `pip install -e ".[docs]"`, GitHub Pages workflow).

### Changed

- **Multiplication kernels no longer expand coefficients.** `grothmult_double` and the
  quantum double kernels keep coefficients as structured products and flat fractions
  instead of calling `expand()`/`cancel()` on every term. This is what makes the quantum
  double Grothendieck kernel usable; `grothmult_double` itself is also considerably
  faster. Call `.simplify()` (or `sympy.cancel`) on the result if a normal form is wanted.
- `BaseRing.from_dict(element)` drops its unused `orig_domain` parameter and builds the
  element in one pass.
- `NilHeckeRing`/`NilHeckeElement` derive from `BaseRing`/`BaseRingElement`.

### Removed

- The research scripts that lived under `src/schubmult/_scripts/` (including the
  `unlinted/` tree and `lr_rule_verify`) are no longer part of the package. Only the eight
  `schubmult_*`/`grothmult_*` CLI entry points remain in `_scripts`.

### Fixed

- `grothmult_double` failed at import on Windows (`resource` is POSIX-only; the memory cap
  is now skipped where it is unavailable).

## 5.0.0

The multiplication kernels are now compiled C++. This release consolidates the
5.0.0b1 and 5.0.0b2 pre-releases.

### Changed

- **The Schubert multiplication kernels run in a C++ extension (`schubmult.schubmult_cpp`).**
  `schubmult_py`, `schubmult_double`, `schubmult_q_fast`, `schubmult_q_double_fast`,
  `schubmult_double_from_elems` and `schubmult_double_alt_from_elems` — and therefore
  every ring product built on them — dispatch to ports of the same algorithms in `cpp/`.
  Results are identical to the Python kernels (coefficients are still unexpanded symengine
  expressions); typical products are one to two orders of magnitude faster. The Python
  kernels remain as the fallback for permutations beyond the compiled size limit
  (`MAXN`, default 32). `SCHUBMULT_NO_CPP=1` forces the Python kernels.
- **Binary wheels for Linux, macOS and Windows** (CPython 3.10–3.14). The extension uses
  only the Python C API — coefficients are handled as `symengine` Python objects — so it
  has no SymEngine C++ dependency and builds from source with any C++17 compiler.
- Ring products no longer re-scan every coefficient's `free_symbols` when assembling the
  result, and products are computed directly in the ring's variables instead of in generic
  variables followed by a substitution pass.
- `requires-python` is now `>=3.10` (3.9 was declared but never supported).
- The PuLP requirement is `PuLP>=2.7.0,<4`. The PuLP 4 pre-releases are a rewrite with an
  incompatible API (no `LpVariable(name=...)`, no bundled CBC) and broke `--display-positive`.

### Added

- Standalone C++ command-line tools in `cpp/` (`schubmult_core`, `schubmult_double_core`,
  `schubmult_q_core`, `schubmult_q_double_core`) mirroring the CLI scripts, including
  `--display-positive` for the double kernel (MILP via `cbc`). These are development tools
  built with `make`; they need the SymEngine C++ headers.

### Fixed

- The script test fixtures are located relative to the test modules, so the test suite
  runs against a non-editable install.

## 4.2.0

### Fixed

- **Parabolic quantum ring constructors were unusable due to a circular import.**
  `schubmult.rings.schubert.schubert_ring` imported `quantum_schubert_ring` at module
  load time, closing a cycle through `parabolic_quantum_double_schubert_ring`. Entering
  the cycle at the parabolic module — as any of `ParabolicQuantumDoubleSchubertRing`,
  `make_parabolic_quantum_basis`, or `make_single_parabolic_quantum_basis` did — raised

  ```
  ImportError: cannot import name 'ParabolicQuantumDoubleSchubertElement' from
  partially initialized module '...parabolic_quantum_double_schubert_ring'
  ```

  The import is now deferred to the single call site that needs it. Reported by an
  external smoke test of 4.1.0, which found the bug while trying to compute QH*(P^3);
  that computation is now a regression test.

- **Five `dual_basis()` declarations raised instead of returning the dual basis.**
  Each imported the right class name from the wrong module:
  `GrothendieckPolyBasis` looked for `GrothendieckBasis` in `free_algebra.schubert_basis`,
  `GlidePolyBasis` and `LascouxPolyBasis` both looked in `free_algebra.fundamental_slide_basis`,
  and `LascouxBasis` looked for `LascouxPolyBasis` in `polynomial_algebra.key_poly_basis`.
  `GrovePolyBasis` declared no dual at all and inherited the base implementation, which
  returns `None`. Every declared dual pair now round-trips in both directions.

### Added

- Grothendieck polynomial multiplication (`schubmult.mult.groth`,
  `schubmult.mult.groth_double`) and a double Grothendieck ring
  (`rings.schubert.double_grothendieck_ring`).
- Lascoux and glide polynomial bases, in both the free algebra
  (`rings.free_algebra.lascoux_basis`, `glide_basis`) and the polynomial algebra
  (`rings.polynomial_algebra.lascoux_poly_basis`, `glide_poly_basis`).
- Chevalley/Monk formula module `rings.schubert.chevalley`.
- Double Grothendieck polynomial separated-descents multiplication (`schubmult.mult.separated_descents`).
- Weak order operations on `Permutation`: `weak_order_leq`, `weak_order_meet`,
  `weak_order_join`.
- Verification scripts `groth_lr_rule` and `verify_quantum_triple_positive`.

### Removed

- `rings.combinatorial.bounded_double_rc_factor_algebra`, an abandoned line of work
  that could not be imported: it required `DecoratedRCGraph`, which had already been
  commented out of `combinatorics.rc_graph`. Its three classes
  (`BoundedDoubleRCFactorAlgebra`, `BoundedDoubleRCFactorAlgebraElement`,
  `BoundedDoubleRCFactorPrintingTerm`) were advertised by `dir(schubmult)` but raised
  `AttributeError` on access. Nothing referenced them. The unrelated
  `BoundedRCFactorAlgebra` is unaffected.
- The commented-out `DecoratedRCGraph` class body and its `full_CEM_double` /
  `elem_squash` helpers, which had no live definitions or callers.

### Notes

- The free algebra is the dual of the polynomial algebra, so each free algebra basis
  is the dual basis to the correspondingly named polynomial basis. The naming
  convention carries the pairing: `SchubertBasis` / `SchubertPolyBasis`,
  `GrothendieckBasis` / `GrothendieckPolyBasis`, `LascouxBasis` / `LascouxPolyBasis`,
  `GlideBasis` / `GlidePolyBasis`, `GroveBasis` / `GrovePolyBasis`, `KeyBasis` /
  `KeyPolyBasis`, `ForestBasis` / `ForestPolyBasis`, `FundamentalSlideBasis` /
  `FundamentalSlidePolyBasis`, `MonomialSlideBasis` / `MonomialSlidePolyBasis`,
  `CompositionSchubertBasis` / `CompositionSchubertPolyBasis`,
  `SeparatedDescentsBasis` / `SepDescPolyBasis`, `ElementaryBasis` /
  `ElemSymPolyBasis`, and `WordBasis` / `MonomialBasis`. Where a dual is itself
  realized inside the free algebra it is named explicitly: `ForestDual`, `GlideDual`,
  `GroveDual`.
- The PyPI project description is generated from `README.md`, which documents the
  current (4.x) API. Releases before this one carried a v1.x-era description.
