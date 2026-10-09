# Changelog

## 5.2.0b2

Pre-release of 5.2.0: SageMath integration (`schubmult.sage`), a correctness fix for parabolic
quantum products with trailing size-1 blocks, and much faster `ElementaryBasis` transitions. Since
5.2.0b1: parabolic quantum Grothendieck products, a manifestly positive expansion of double
Grothendieck polynomials in double Schubert polynomials, a transition + Monk kernel (lrcalc's
algorithm) with a cost-model hybrid, probabilistic zero testing for the double and quantum double
kernels, and fast normalization of double Grothendieck coefficients. Pre-releases install with
`pip install --pre schubmult`; `pip install schubmult` keeps giving 5.1.1 until 5.2.0 is final.

### Added

- **Positive expansion of double Grothendieck polynomials in double Schubert polynomials**
  (`schubmult.mult.groth_double.dgroth_to_dschub_positive`). Writes
  $\mathfrak G_v(x;z)=\sum_w a_{v,w}(\beta,z)\,\mathfrak S_w(x;z)$ with every $a_{v,w}$ in
  $\mathbb N[\beta,z]$, as a sum over triples of pipe dreams with no sign and no change of variables:
  the Fomin--Kirillov Cauchy factorization over Demazure products
  $\mathfrak G_v(x;z)=\sum_{u\,\star\,v'=v}\beta^{\ell(u)+\ell(v')-\ell(v)}\mathfrak G_{u^{-1}}(z)\mathfrak G_{v'}(x)$,
  then Lenart's integer expansion of $\mathfrak G_{v'}(x)$ in Schubert polynomials
  (`WCGraph.groth_to_schub`), then $\mathfrak S_{u'}(x)=\sum_{u'=aw}\mathfrak S_a(z)\mathfrak S_w(x;z)$.
  Agrees with the exact divided-difference expansion `dgroth_to_dschub` on all of $S_4$ and samples
  in $S_5$. Also new, and conjectural: `dgroth_positive_phantom_expansion` (the sign-normalized
  staircase phantom expansion in the formal-inverse alphabet, with positivity checked at every step
  and a `ValueError` on failure) and `dgroth_copipe_expansion`, an unsigned reduced-pipe-dream /
  co-pipe-dream formula for the same coefficients, verified against the recurrence for every
  permutation in $S_6$. These are small-rank research routines, not multiplication kernels.
- **Parabolic quantum Grothendieck products.** `grothmult_q` and `grothmult_q_double` accept
  `--parabolic g1 g2 ...` (block sizes) and return the product in $QK(G/P)$ resp. $QK_T(G/P)$ for the
  partial flag variety $Fl(g_1, g_1+g_2, \ldots; n)$; inputs must be minimal coset representatives.
  Uses Kato's ring homomorphism $QK_T(G/B)\twoheadrightarrow QK_T(G/P)$,
  $\mathcal O^w\mapsto\mathcal O^{[w]_P}$, $Q_j\mapsto 1$ for $j\in P$ (arXiv:1906.09343, Thm 2.19),
  so unlike `schubmult_q --parabolic` there is no Peterson--Woodward comparison term by term. Exposed
  on the web app too.
- **Transition + Monk kernel and a cost-model hybrid for single Schubert products**
  (`schubmult.mult.transition`). `schubmult_py_transition` multiplies by expanding one factor into
  monomials with the Lascoux--Schützenberger transition recursion (`transition_monomials`) and
  folding them into the other factor one variable at a time by Monk's rule, as a Horner scheme over
  the Schubert basis -- lrcalc's algorithm. Its cost is driven by the number of pipe dreams
  (`pipe_dream_count`), the v-path kernel's by the number of v-paths, and an exhaustive comparison
  over $S_7$ (all 12.7M unordered pairs) and scaling families to $S_{12}$ showed each is faster by
  large factors where the other is slow; `schubmult_py_hybrid` picks between them with a cost model
  fitted on that data. Both have C++ implementations in the extension. `schubmult_py` and `Sx` are
  unchanged; the hybrid is opt-in.
- **Probabilistic zero testing in the double and quantum double kernels**
  (`schubmult_double --probabilistic`, `schubmult_q_double --probabilistic`;
  `DoubleSchubertPolynomialRing(R, probabilistic=True)` in Sage). The layered DP of
  `schubmult_double` has signed terms that cancel: a state whose partial sum is zero as a polynomial
  but not structurally keeps fanning out, and its leaves are output terms with coefficient zero
  (4187 terms returned for a 246-term answer). Partial sums are now shadowed by exact evaluations
  at random integer points (Schwartz--Zippel; `schubmult.mult._shadow`) and dead states are pruned
  as they arise. `X([4,1,6,5,2,3]) * X([8,1,7,6,2,3,5,4])` in Sage went from 11 s to 0.5 s. The
  quantum double kernel keys its states by $q$-monomial, which alone makes the exact path ~18%
  faster; with `--probabilistic` the benchmark `--code 2 0 5 5 0 5 - 3 0 5 0 3 3` halves again
  (33 s / 2.65 GB -> 14.5 s / 0.84 GB). Exposed on the web app as a checkbox.
- **Unexpanded coefficients in the Sage rings** (`raw_coefficients=True` on
  `DoubleSchubertPolynomialRing`, `GrothendieckPolynomialRing`, ...). Every native Sage polynomial
  ring stores expanded normal forms, and converting the kernels' factored coefficients through
  libsingular dominated the time of large products. With `raw_coefficients=True` the base ring is
  `schubmult.sage.symengine_ring.SymEngineRing`, whose elements are the kernel's SymEngine
  expressions as is; zero testing is probabilistic (`nonzero_mask`, `is_identically_zero`) and
  printing normalizes all coefficients at once with shared subtrees.
- **`--simplify` on `grothmult_double`, `grothmult_q_double`** (and `schubmult_double`,
  `schubmult_q_double`): normalize the output coefficients. Slower, but the printed result shrinks
  from tens of thousands of characters per coefficient to something readable.
- **Web app: download results as a file.** A checkbox emits the full result as a file download
  instead of the text box, with its own limit (`SCHUBMULT_MAX_DOWNLOAD_BYTES`, default 100 MB)
  separate from the display limit (`SCHUBMULT_MAX_OUTPUT_BYTES`, default 1 MB), for users who
  actually want a large product rather than a truncated one.
- **Baseline type hints** (`schubmult._typing`: `PermLike`, `PermCoeffDict`, `Alphabet`, `Expr`,
  `Coeff`) on `Permutation`, `RCGraph`, `WCGraph`, the Grothendieck kernels and the variable
  machinery, checked with mypy on that subset.
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
  - `GrothendieckPolynomialRing(R)` and `DoubleGrothendieckPolynomialRing(R, alphabet='y')` --
    $\beta$-Grothendieck polynomials $\mathfrak G^\beta_w(x)$ over `R[beta]` and their double
    versions $\mathfrak G^\beta_w(x; y)$ (with $x \oplus y = x + y + \beta x y$) over the fraction
    field of `R[beta, y, z]`, where the equivariant K-theoretic structure constants live (denominators
    are products of $1 + \beta y_i$). `beta()` returns the deformation parameter; $\beta = 0$
    recovers the (double) Schubert polynomials. Coefficients of the double ring are reduced
    rational functions computed without any symbolic expansion on the schubmult side (a flat
    numerator / atom-power denominator evaluation in the libsingular ring underneath, one gcd per
    coefficient); a triple product of `S_4` elements with 24 terms and numerators of ~10^4 terms
    takes a few seconds.
  - `PolynomialAlgebra(R)` and `PolynomialAlgebra(R, n)` -- the polynomial ring `R[x0, x1, ...]`
    (resp. in `n` variables) as a Sage parent with realizations (like `SymmetricFunctions`), each
    indexed as in the literature: `schubert()` and `grothendieck()` by permutations (`S[3, 1, 2]`);
    `monomial()`, `key()`, `fundamental_slide()`, `monomial_slide()`, `forest()` and the K-theoretic
    `glide()`, `lascoux()`, `grove()` by weak compositions (`k[2, 0, 1]`); and, in `n` variables only,
    `elementary()` -- products of elementary symmetric polynomials `e_a(x0..x_{j-1})`, whose indexing
    depends on `n`. The K-theoretic bases are at `beta = -1`, the classical convention
    (`G[1, 3, 2].expand()` is `x0 + x1 - x0*x1`); `GrothendieckPolynomialRing` and
    `DoubleGrothendieckPolynomialRing` are unaffected and keep `beta` as a variable of their base ring.
    schubmult's own bases work at `beta = 1` and the
    layer converts through the grading (`P^{-1}_a(x) = (-1)^|a| P^{1}_a(-x)`). All bases coerce into one another through the
    monomial basis, products come from the schubmult bases and Schubert kernels, Sage's
    `SchubertPolynomialRing` and `KeyPolynomials` coerce in, and the elements coerce into the (double,
    quantum, Grothendieck) Schubert rings above. schubmult's own `PolynomialAlgebra` is graded by the
    number of variables (keys of different lengths multiply to zero, the structure dual to the free
    algebra); the Sage parents are the plain polynomial ring and pad keys to a common number of
    variables before calling schubmult.
  - All accept lists/`Permutation`s, Sage polynomials (finite or `InfinitePolynomialRing`), and
    elements of `SchubertPolynomialRing`, `KeyPolynomials` and `AtomPolynomials` (whose variables are
    read as `x`, as Sage's Schubert ring does; these also coerce, so `X(w) + k([1])` works).
    `project(n)` drops the classes that vanish in the flag variety of `C^n` (the rings are stable, so
    a product carries every class of the infinite flag variety; e.g. Buch's equivariant quantum
    `Gr(2, 5)` example is `(QD(u) * QD(v)).project(5)` with `parabolic=(2, 3)`).
    `from_symmetric_function(f, n)` expands a symmetric function in `x0..x{n-1}` (a Schur function
    gives the Grassmannian Schubert polynomial with descent at `n`) and `to_symmetric_function(n=None)`
    reads a symmetric element back as a Schur expansion over the coefficient variables (Grassmannian
    double Schubert polynomials are factorial Schur functions). `expand()` lands in
    `R[x0.., y0.., q0..]`, 0-indexed like Sage; `divided_difference(i)` on the double ring.
    Coefficients are converted at the boundary, so nothing SymEngine-flavoured leaks into the Sage API.
  - Sage is **not** a dependency: the subpackage is inert unless imported, and
    `tests/test_api_surface.py::test_package_is_fully_usable_without_sage` runs the whole package
    with Sage imports blocked. Install schubmult into a Sage environment
    (`sage -pip install schubmult`) to use it; CI runs the Sage doctests in a separate job against
    conda-forge Sage.
- **`strip_zeros(exact=True)` on every ring element.** The kernels leave coefficients as products of
  factors, and in large products a sizeable fraction of them cancel to zero without looking like it
  (7868 -> 3845 terms for a quantum double product in $S_8$). `strip_zeros()` (now on
  `BaseRingElement`, so available in every ring) still only drops literal zeros; `exact=True` also drops
  coefficients that are zero as polynomials, detected by evaluating them at random integer points in
  exact arithmetic with memoized shared subtrees (`schubmult.symbolic.functions.vanish_at_random_points`)
  -- about 0.7x the cost of the product itself, versus ~10x for expanding the coefficients. Products
  are unchanged: stripping stays opt-in.

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
- **`CompositionSchubertPolyBasis` received `(perm, length)` keys.** Transitions landing in a
  Schubert-type basis (from `MonomialBasis`, `GrothendieckPolyBasis`, and `SchubertPolyBasis` itself)
  handed the target the keys of `SchubertPolyBasis` without normalizing them, so
  `PolynomialAlgebra(CompositionSchubertPolyBasis(x)).from_expr(...)` and every `change_basis` into
  that basis produced elements keyed by permutations instead of Lehmer codes, and any transition out
  of the ring then failed. Keys now go through the target basis's `attach_key`.

### Changed

- **Double Grothendieck coefficients are normalized in the multiplicative variables.**
  `DoubleGrothendieckElement.simplify()` ran `sympy.cancel` on the raw kernel output; coefficients
  of one product in $K_T(Fl_5)$ reach ~65k characters unexpanded and took 8--58 s each (274 s for
  `DGx([5,1,3,2,4]) * DGx([2,1,5,4,3])`, 56 terms; a 10-product run was 332 s against 0.2 s for
  Buch's EquivCalc). Cancellation is now done in $t_i = 1+\beta y_i$: every coefficient is a Laurent
  polynomial in the $t_i$, the denominator atoms become monomials, a single expansion collects
  everything and the denominator is read off as the monomial of negative exponents. No gcd is
  computed and nothing is expanded in $y$.
- **`grothmult_double` keeps its result factored as Grothendieck factorial elementary symmetric
  functions.** The kernel now parallels `schubmult_double` directly -- the v-path's layer factors are
  the $K$-theoretic factorial elementaries `groth_elem_sym_poly` / `groth_elem_sym_func`, with the
  Pieri coefficients in closed form (`_tilde_elem_sym_frac`) -- instead of transitioning
  $\mathfrak G_v$ to the double Schubert basis first and applying the Schubert Pieri rule. The
  quantum double kernel iterates its v-paths top-down for the same reason. `BoundedWCFactorAlgebra`'s
  `full_groth_elem` likewise runs the Grothendieck v-path directly, with each layer factor
  $\tilde E_{p,k}(x;0)=G_{1^p}(x_1..x_k)+\beta\,G_{1^{p+1}}(x_1..x_k)$ as a sum of elementary
  Grothendieck WC graphs (`groth_elem_factor`), rather than approximating with Schubert elementaries.
- **`positivity.py`** lost years-old, misleading comments and picked up a few small optimizations;
  behaviour of `--display-positive` is unchanged.
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
  .[test]` locally) and triggers on pull requests into `develop` as well as `main`.
- **Integration branch renamed `redevelop` -> `develop`.** Feature branches PR into `develop`;
  `main` only receives final releases from `develop` or a `hotfix/*` branch.
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
