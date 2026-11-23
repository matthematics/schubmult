<a id="schubmult.mult.positivity"></a>

# schubmult.mult.positivity

Manifestly positive representations of double Schubert and Grothendieck structure constants.

The structure constants ``c^w_{u,v}(y, z)`` of double Schubert polynomial
multiplication (``schubmult_double``) are known to be polynomials in the
differences ``y_i - z_j`` with nonnegative integer coefficients (Graham's
positivity theorem). This module computes that manifestly positive form:

- ``posify``: the main recursive engine. Reduces ``(u, v, w)`` via known
  combinatorial identities (pattern-avoidance checks, dominance/one-dominance,
  descent and coefficient reductions in ``schubmult.utils.schub_lib``) down to
  cases with closed positive formulas (``dualcoeff``, ``forwardcoeff``, or a
  single elementary symmetric polynomial), falling back to the integer-LP
  solver ``compute_positive_rep`` when no reduction applies.
- ``compute_positive_rep``: expresses an arbitrary such polynomial as a
  nonnegative-integer combination of product-of-differences monomials, found
  via an integer program (PuLP) over a spanning set of candidate monomials.
- ``dualcoeff``/``forwardcoeff``/``dualpieri``: closed-form positive rules for
  special cases (``u`` dominates ``w``, the ``will_formula_work`` forward Monk
  case, and the dual Pieri expansion respectively).
- ``groth_posify``: positive FGL representations of double Grothendieck
  coefficients, using an integer program and restoring the symbolic ``beta`` grading.
  ``groth_posify_many`` does a whole batch of coefficients, solving their integer
  programs concurrently.

<a id="schubmult.mult.positivity.cbc_solver"></a>

#### cbc\_solver

```python
def cbc_solver(msg=False)
```

PuLP ``COIN_CMD`` using the CBC binary bundled by ``pulp[cbc]``, falling back to ``cbc`` on PATH.

<a id="schubmult.mult.positivity.groth_posify"></a>

#### groth\_posify

```python
def groth_posify(val, var2, var3, msg=False)
```

Positive FGL form of a coefficient, at ``beta = 1``.

Basis: all products of ``z_i - y_j``, ``(1 + y_j)^{-1}``, ``(1 + z_i)^{-1}``.
Substituting ``u_j = 1 + y_j``, ``v_i = 1 + z_i`` (so ``y = u - 1``,
``z = v - 1``) turns the value into a Laurent polynomial in ``u, v``; clearing
the monomial denominator makes everything polynomial.  Candidates are
difference products times leftover clearing monomials, expanded once, with
``as_coefficients_dict`` terms used as opaque basis vectors for an integer
LP.  ``beta`` is restored as ``beta**(k - d)`` for a term with ``k`` difference
factors, where ``d`` is the homogeneous degree of the value (``deg y = deg z = 1``,
``deg beta = -1``), the Laurent atoms ``(1 + beta*y)`` having degree 0 by construction.

The problem splits into independent blocks by total degree.  Each block is
reduced by an exact presolve and the remainder is solved
as an integer program, escalating the number of difference factors.  The
result is verified exactly: by an integer identity at ``beta = 1`` when
homogeneity is structurally evident, otherwise by symbolic expansion.

``msg`` is passed through to the solver as its verbosity option.  Use
:func:`groth_posify_many` for several coefficients, whose integer programs can
then be solved concurrently.

<a id="schubmult.mult.positivity.groth_posify_many"></a>

#### groth\_posify\_many

```python
def groth_posify_many(vals,
                      var2,
                      var3,
                      msg=False,
                      *,
                      max_workers=None,
                      progress=None)
```

:func:`groth_posify` of each coefficient in ``vals``, solving the integer programs concurrently.

The solver subprocesses are driven by up to ``max_workers`` threads, while the symbolic work and the construction
caches shared between coefficients stay on the calling thread.  The results do not depend on ``max_workers``: they
are those of sequential :func:`groth_posify` calls, yielded in the order of ``vals`` as soon as each is ready.  If a
coefficient has no positive representation, the results before it are yielded and then its exception is raised.
Closing the generator early, or interrupting it, stops the solver processes and removes their temporary files.

**Arguments**:

- `vals` - Coefficients, as for :func:`groth_posify`.
- `var2` - First secondary alphabet (``y``).
- `var3` - Second secondary alphabet (``z``).
- `msg` - Passed through to the solver as its verbosity option.
- `max_workers` - Number of integer programs solved at once, at least 1; 1 runs everything on the calling
  thread.  Defaults to 1.  If ``msg`` is enabled with multiple workers, solver logs may interleave.
  Each solver process needs memory of its own, which is large for the hardest cases.
- `progress` - Called as ``progress(i)``, on the calling thread, when work on ``vals[i]`` starts.
  

**Returns**:

  A generator of the ``groth_posify`` results, in the order of ``vals``.
  

**Raises**:

- `ValueError` - If ``max_workers`` is less than 1.

<a id="schubmult.mult.positivity.compute_positive_rep"></a>

#### compute\_positive\_rep

```python
def compute_positive_rep(val, var2=None, var3=None, msg=False)
```

Express ``val`` as a nonnegative-integer combination of product-of-differences monomials.

``val`` must be a polynomial in ``var2``/``var3`` known (by positivity of
double Schubert structure constants) to admit an expansion
``sum_b n_b * prod (var2_i - var3_j)`` with ``n_b >= 0`` integers. Builds a
candidate spanning set of such product monomials from ``val``'s own
monomials, then solves an integer program (via PuLP) for nonnegative
integer coefficients ``n_b`` matching ``val`` exactly.

**Arguments**:

- `val` - Symbolic polynomial expression in ``var2``/``var3``.
- `var2` - First secondary alphabet (``y``).
- `var3` - Second secondary alphabet (``z``).
- `msg` - Passed through to the LP solver as its ``msg`` (verbosity) option.
  

**Returns**:

  A symbolic expression equal to ``val``, written as a sum of
  nonnegative-integer multiples of product-of-differences monomials.
  

**Raises**:

- `Exception` - If the reconstructed expression does not equal ``val``
  (i.e. no valid nonnegative integer solution reproduces it exactly).

<a id="schubmult.mult.positivity.posify"></a>

#### posify

```python
@cached(
    cache={},
    key=lambda val, u2, v2, w2, var2=None, var3=
    None, msg=False, optimize=True, n=_vars.n: hashkey(
        val, u2, v2, w2, var2, var3, msg, optimize, n),
)
def posify(val,
           u2,
           v2,
           w2,
           var2=None,
           var3=None,
           msg=False,
           optimize=True,
           n=_vars.n)
```

Manifestly positive representation of the structure constant ``c^{w2}_{u2,v2}(var2, var3)``.

``val`` is the (already computed, possibly not manifestly positive) value of
the coefficient of ``S_{w2}`` in ``S_{u2}(x, var2) * S_{v2}(x, var3)``.
Recursively reduces ``(u2, v2, w2)`` via pattern-avoidance-guarded identities
(``try_reduce_u``/``try_reduce_v``, ``reduce_descents``, ``reduce_coeff``,
``is_split_two``) toward cases handled by closed positive formulas
(a single elementary symmetric polynomial when ``v`` has one nonzero code
entry, ``dualcoeff`` when ``will_formula_work(v, u)`` or ``u`` dominates
``w``, ``forwardcoeff`` when ``will_formula_work(u, v)``, or the
length-one-difference case built from ``pull_out_var``/``schubpoly``
directly). Falls back to ``compute_positive_rep`` (an integer-LP search)
when no reduction or closed formula applies and ``optimize`` is true.

Results are cached by all arguments.

**Arguments**:

- `val` - The structure constant to re-express positively.
- `u2` - First factor's permutation.
- `v2` - Second factor's permutation.
- `w2` - Target permutation (coefficient of ``S_{w2}``).
- `var2` - First secondary alphabet.
- `var3` - Second secondary alphabet.
- `msg` - Verbosity flag passed down to ``compute_positive_rep``'s LP solver.
- `optimize` - If ``False``, return ``val`` unchanged when no closed-form
  reduction applies (skip the LP fallback); if ``None`` and that
  case is reached, raise.
- `n` - Size of the ambient alphabet used when no other bound is available.
  

**Returns**:

  A manifestly positive expression equal to ``val``, or ``val`` itself
  when no formula applies and optimization is disabled.

<a id="schubmult.mult.positivity.shiftsub"></a>

#### shiftsub

```python
def shiftsub(pol, var2=None)
```

Shift every ``var2[i]`` in ``pol`` up to ``var2[i + 1]`` (for ``i`` in ``0..98``).

<a id="schubmult.mult.positivity.posify_generic_partial"></a>

#### posify\_generic\_partial

```python
def posify_generic_partial(val, u2, v2, w2)
```

``posify`` specialized to the fixed generic alphabets ``_vars.var_g1``/``_vars.var_g2``.

Asserts (raises on mismatch) that the recomputed positive expression equals
the input ``val``, as a consistency check.

<a id="schubmult.mult.positivity.schubmult_generic_partial_posify"></a>

#### schubmult\_generic\_partial\_posify

```python
@cache
def schubmult_generic_partial_posify(u2, v2)
```

Manifestly positive expansion of ``S_{u2}(x, var_g1) * S_{v2}(x, var_g2)``.

Returns ``{w2: coeff}`` where each ``coeff`` is the positive representation
(via ``posify_generic_partial``) of the corresponding
``schubmult_double_pair_generic_alt`` coefficient.

<a id="schubmult.mult.positivity.forwardcoeff"></a>

#### forwardcoeff

```python
def forwardcoeff(u, v, perm, var2=None, var3=None)
```

Closed-form structure constant ``c^{perm}_{u,v}(var2, var3)`` for the "forward" case
(used when ``will_formula_work(u, v)`` holds in ``posify``).

Writes ``muv = uncode(v.theta())`` and reduces to a lookup in
``schubmult_double_pair(u, muv, var2, var3)`` when the length condition
``(perm * (~v * muv)).inv == (~v * muv).inv + perm.inv`` holds; returns 0 otherwise.

<a id="schubmult.mult.positivity.dualcoeff"></a>

#### dualcoeff

```python
def dualcoeff(u, v, perm, var2=None, var3=None)
```

Closed-form structure constant ``c^{perm}_{u,v}(var2, var3)`` for the "dual" case
(used in ``posify`` when ``will_formula_work(v, u)`` holds or ``u`` dominates ``perm``).

When ``u`` is the identity, reduces directly to a single Schubert
polynomial ``schubpoly(v * (~perm), var2, var3)``. Otherwise expands via
``dualpieri`` (directly if ``u`` dominates ``perm``, or after rewriting to
``u``'s dominant permutation ``uncode(u.theta())`` otherwise), summing
products of ``(var2_{i+1} - var3_j)`` factors against a final
``schubpoly`` term.

<a id="schubmult.mult.positivity.dualpieri"></a>

#### dualpieri

```python
def dualpieri(mu, v, w)
```

Dual Pieri expansion used by ``dualcoeff``: enumerate the data witnessing
``S_mu * S_v -> S_w`` when ``mu`` is dominant.

Compares ``mu``'s inverse code against ``w``'s inverse code layer by layer,
peeling one "cycle" of variables per layer via ``divdiffable``/``pull_out_var``,
and returns the list of ``[vlist, vp]`` pairs consumed by ``dualcoeff`` to
build the final positive expression (empty list if ``w`` is not reachable
from ``mu``, ``v`` this way).

