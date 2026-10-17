<a id="schubmult.mult.positivity"></a>

# schubmult.mult.positivity

Manifestly positive representations of double Schubert structure constants.

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

<a id="schubmult.mult.positivity.cbc_solver"></a>

#### cbc\_solver

```python
def cbc_solver(msg=False)
```

PuLP ``COIN_CMD`` using the CBC binary bundled by ``pulp[cbc]``, falling back to ``cbc`` on PATH.

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
    None, msg=False, sign_only=False, optimize=True: hashkey(
        val, u2, v2, w2, var2, var3, msg, sign_only, optimize),
)
def posify(val,
           u2,
           v2,
           w2,
           var2=None,
           var3=None,
           msg=False,
           sign_only=False,
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

Results are cached by ``(val, u2, v2, w2, var2, var3, msg, sign_only, optimize)``.

**Arguments**:

- `val` - The structure constant to re-express positively.
- `u2` - First factor's permutation.
- `v2` - Second factor's permutation.
- `w2` - Target permutation (coefficient of ``S_{w2}``).
- `var2` - First secondary alphabet.
- `var3` - Second secondary alphabet.
- `msg` - Verbosity flag passed down to ``compute_positive_rep``'s LP solver.
- `sign_only` - If ``True``, only determine and return the sign of ``val``
  (``-1``, ``0``, or ``1``) rather than a full positive expression.
- `optimize` - If ``False``, return ``val`` unchanged when no closed-form
  reduction applies (skip the LP fallback); if ``None`` and that
  case is reached, raise.
- `n` - Size of the ambient alphabet used when no other bound is available.
  

**Returns**:

  A manifestly positive expression equal to ``val`` (or, if
  ``sign_only``, one of ``-1``, ``0``, ``1``).

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

