# Changelog

## 5.2.0 (in pre-release: 5.2.0b1)

SageMath integration (`schubmult.sage`), a correctness fix for parabolic quantum products with
trailing size-1 blocks, and much faster `ElementaryBasis` transitions. Pre-releases install with
`pip install --pre schubmult`; `pip install schubmult` keeps giving 5.1.1 until 5.2.0 is final.

### Added

- **`schubmult.sage`: SageMath integration.** Sage parents built on `CombinatorialFreeModule`
  whose arithmetic is delegated to the schubmult kernels, alongside Sage's own
  `SchubertPolynomialRing`:
  - `DoubleSchubertPolynomialRing(R, alphabet='y')` -- double Schubert polynomials
    $\mathfrak S_w(x; y)$ over `R[y, z]`; a ring in another alphabet
    (`DoubleSchubertPolynomialRing(R, 'z')`) coerces in, so `X(u) * Z(v)` is the mixed product
    $\mathfrak S_u(x;y)\,\mathfrak S_v(x;z)$ expanded in the `y`-basis.
  - `QuantumSchubertPolynomialRing(R, parabolic=None)` and
    `QuantumDoubleSchubertPolynomialRing(R, alphabet='y', parabolic=None)` -- quantum and
    quantum double Schubert polynomials over `R[q]` resp. `R[q, y, z]`; `parabolic=(n_1, ..., n_k)`
    gives the partial flag variety with those block sizes, with an implicit unbounded last block
    (basis permutations may have descents exactly at the recorded boundaries and be increasing
    beyond them; products extend by one new block as needed, never by enlarging a recorded one).
  - All accept lists/`Permutation`s, Sage polynomials (finite or `InfinitePolynomialRing`), and
    elements of `SchubertPolynomialRing`; `expand()` lands in `R[x0.., y0.., q0..]`, 0-indexed
    like Sage; `divided_difference(i)` on the double ring. Coefficients are converted at the
    boundary, so nothing SymEngine-flavoured leaks into the Sage API.
  - Sage is **not** a dependency: the subpackage is inert unless imported, and
    `tests/test_api_surface.py::test_package_is_fully_usable_without_sage` runs the whole package
    with Sage imports blocked. Install schubmult into a Sage environment
    (`sage -pip install schubmult`) to use it; CI runs the Sage doctests in a separate job against
    conda-forge Sage.

### Fixed

- **Parabolic quantum products with a trailing size-1 block** (`QPSx(2, 3, 1)`,
  `QPSx(2, 3, 1, 1)`, ...). `apply_peterson_woodward` inferred the ambient flag size from the last
  parabolic generator, and a size-1 block contributes none, so the ambient was taken one block
  short: result permutations of length exactly `sum(blocks)` were dropped (for
  $\sigma_{32}\sigma_1$ with blocks `(2,3,1)` the class $\sigma_{42}$ vanished) and the product no
  longer matched the basis polynomials. `apply_peterson_woodward(..., n=)` now takes the ambient
  size explicitly and the parabolic rings pass `sum(blocks)`. `QPSx(2, 3)`, `(2, 3, 1)`,
  `(2, 3, 2)` and `(2, 3, 1, 1)` now all agree, as they should. The CLI (`--parabolic`), which only
  knows generator indices, is unchanged.

### Changed

- **`ElementaryBasis` <-> `SchubertBasis` transitions are computed blockwise.** Both directions
  go through the finite `(numvars, degree)` block: every elementary product of that degree is
  expanded in the Schubert basis by the Pieri rule (`ElementaryBasis.schubert_block`, cached), which
  is the Schubert -> Elem matrix, and its inverse is Elem -> Schubert. This replaces
  `Sx.from_expr` on staircase-padded monomial-symmetric polynomials and RC-graph enumeration of
  `S_{perm w0}`; e.g. `FA(2) * FA(2, 1)` in the elementary basis went from 25 s to 0.04 s. Results
  are identical on every key checked ($n \le 3$, degree $\le 4$). `ElementaryBasis.degree_keys`
  and `transition_from_schubert` are new; `ElementaryBasis.staircase` is gone.
- **Test suite runs in ~1 min instead of ~8.** Slow tests were made cheaper without losing
  coverage: exact zero tests of rational-function identities now evaluate at random rational
  points (`schubmult.utils.test_utils.vanishes`) instead of `sympy.cancel` on huge unsimplified
  differences (90 s -> 0.2 s for the double Grothendieck round trip); the heaviest
  `--display-positive` script fixtures were replaced by same-flag, smaller permutations; a few
  ring examples dropped one degree. CI runs `pytest -n auto` (`pytest-xdist`; `pip install -e
  .[test]` locally) and triggers on pull requests into `redevelop` as well as `main`.
- **Docs deploy only from final release tags.** The docs workflow no longer runs on pushes to
  `main`; it runs on release tags and refuses `.dev`, local (`+...`), and pre-release
  (`a`/`b`/`rc`) versions, so a pre-release tag publishes wheels but leaves the documentation site
  at the last final release.

## 5.1.1

Patch release: compatibility with PuLP 4.0.0, a correctness fix for the elementary basis of
the free algebra, and a much lighter import footprint (SymPy loads only when needed).

### Fixed

- **`--display-positive` under PuLP 4.0.0.** PuLP 4.0.0 (released after 5.1.0) broke
  `compute_positive_rep` and `grothmult_double --display-positive` in several ways:
  `expr == <symengine Integer>` now yields a bare `False` instead of a constraint
  (`TypeError: A False object cannot be passed as a constraint`), `pulp.LpStatus` was removed
  and `LpProblem.solve` returns an `LpSolveStats` object, `LpAffineExpression(dict)` can no
  longer form constraints, and variables become unusable once their `LpProblem` is garbage
  collected. Constraints are now built with `int(...)` right-hand sides and `lpSum`, status is
  read via a helper that accepts both APIs, and variable values are read while the problem is
  alive. The code runs on both PuLP 3.3.x and 4.0; the requirement is
  `PuLP[cbc]>=3.3.2, <4` for now so that a resolver does not pull in 4.x untested.
- **`ElementaryBasis` / `ElemSymPolyBasis` were not dual.** Keys are `(tup, numvars)` with
  `tup[:numvars-1]` the flag part ($e_{a_i}(x_1..x_i)$) and `tup[numvars-1:]` the symmetric
  tail (a sorted product of $e_k(x_1..x_n)$). `ElementaryBasis.transition_schubert` built its
  staircase from the key's own tail length, so `Elem(key)` only annihilated the `E(key')`
  whose tail fit inside that staircase; e.g. $\langle \mathrm{Elem}((0,2),2),\,
  E((0,1,1),2)\rangle = 1$ ($e_2$ against $e_1^2$). The staircase is now padded uniformly to
  the degree, `SchubertBasis.transition_elementary` uses the same staircase and emits
  canonical keys (zeros stripped from the tail, tail sorted, `(0,)` if empty), and the
  `numvars == 1` slicing bug in `transition_schubert` is gone. `Elem -> Schub -> Elem` is
  now the identity and the delta property holds on every canonical key checked
  ($n \le 4$, degree $\le 4$). `SchubertPolyBasis.transition_elementary` emits the same
  canonical keys.
- `MonomialSlideBasis` / `MonomialSlidePolyBasis` declare each other as `dual_basis()`;
  `GrovePolyBasis` expands through its monomial basis so `expand()` returns a polynomial;
  `GrothendieckPolyBasis` products go through the `grothmult_py` kernel at $\beta = 1$
  instead of the ring-level route. The polynomial-algebra coproduct is the free-algebra
  product's adjoint again for every declared dual pair.
- `BaseSchubertRing` defines `__hash__` consistent with its `__eq__` (type and generating
  sets), so rings can be dict keys and cache keys; ring *elements* are explicitly unhashable
  (`__hash__ = None`) since they are mutable dicts.

### Changed

- **SymPy is imported lazily.** `schubmult.symbolic` exposes SymPy names as
  `schubmult.utils._lazy.LazyAttr` proxies that resolve on first call, attribute access or
  subclassing; the ring-domain protocol (`EXRAW`, `CoercionFailed`, `Ring`, ...) lives in the
  new SymPy-free `schubmult.symbolic.domain`; ring elements print through
  `schubmult.utils._printable.LazyPrintable` instead of subclassing SymPy's `Printable`;
  `schubmult.mult` and the kernels in `schubmult.mult.single` import SymPy/SymEngine only
  inside the functions that need them. `schubmult_py` with integer output never loads SymPy,
  and the CLI sets `OPENBLAS_NUM_THREADS`/`OMP_NUM_THREADS`/`MKL_NUM_THREADS=1` before
  importing numpy. Import time and CLI start-up drop accordingly; no public names moved.
- Ring elements print as linear text built from SymEngine's `str` (terms sorted by length
  then permutation) instead of SymPy's pretty printer, which took seconds on large double
  Grothendieck coefficients. `_repr_latex_` is disabled for these elements so notebooks show
  the fast form; call `latex(elem)` or `pretty(elem)` explicitly for the SymPy renderings.
- The PuLP requirement is `PuLP[cbc]>=3.3.2, <4` (was `PuLP>=2.7.0,<4`). `--display-positive`
  builds variables with `LpProblem.add_variable` and solves with `COIN_CMD` using the CBC
  binary from the `[cbc]` extra (`cbcbox`), replacing `LpVariable(...)` and `PULP_CBC_CMD`.
- Ruff `UP038` is no longer ignored; `isinstance` checks use `X | Y` unions.

### Added

- `benchmark_schubmult.py`: timing harness for the `schubmult_*`/`grothmult_*` kernels.
- `web/`: a small Flask wrapper exposing the CLI scripts through a web form and
  `POST /api/compute`, embeddable via `<iframe>` (see `web/README.md`, `web/DEPLOY.md`).
- Tests: `tests/rings/test_duality.py` covers `ElementaryBasis` and `MonomialSlideBasis`,
  includes colliding elementary keys, and checks that the free-algebra product is adjoint to
  the polynomial coproduct (and vice versa) for every dual pair; more polynomial-algebra and
  free-algebra basis round-trip tests in `tests/rings/`.
- CI runs on pushes to `redevelop` as well as `main`.

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
