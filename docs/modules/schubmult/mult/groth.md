<a id="schubmult.mult.groth"></a>

# schubmult.mult.groth

Multiplication kernels for (single) beta-Grothendieck polynomials.

* ``grothmult_py``: the product ``G_u * G_v`` in the ``G`` basis, by the same route as
  ``groth_double.grothmult_double`` specialized to ``y = z = 0`` -- expand ``G_v`` into
  Schubert polynomials via WC graphs, then push each ``S_{v'}`` through the v-path layers
  of ``theta(v'^{-1})`` with the closed-form K-Pieri coefficient ``groth_elem_sym_coeff``
  (a binomial times a power of ``beta``).
* ``single_variable_groth`` / ``mult_poly_groth``: multiplication by ``x_k`` (the
  non-equivariant K-theoretic Chevalley formula) and by arbitrary polynomials in ``x``.

The Chevalley coefficients are derived from M. Willems, "A Chevalley formula in
equivariant K-theory" (arXiv:math/0603220), Theorem 5 (the ordinary, non-equivariant
specialization of his equivariant Chevalley formula, Theorem 4). Willems indexes K-theory
classes O_w by the *dimension* of the Schubert variety, dual to the *codimension* indexing
used by Schubert/Grothendieck polynomials S_w/G_w; the w0-conjugation in
``_chevalley_ev_weights`` (``hat_w = w0*w`` going in, ``w0*v`` coming out) translates
between the two conventions. The beta-grading (beta^(d-1) per length difference
d = l(v)-l(w)) matches this codebase's beta-deformed Grothendieck polynomial
normalization (beta=0 recovers the classical Monk formula).

<a id="schubmult.mult.groth.chevalley_x_k"></a>

#### chevalley\_x\_k

```python
def chevalley_x_k(w, k, beta, n=None)
```

Coefficients of ``x_k * G_w^(beta)`` in the Grothendieck basis, as a dict
``{v: coeff}`` (``w`` itself never appears: the self-term cancels identically).

<a id="schubmult.mult.groth.single_variable_groth"></a>

#### single\_variable\_groth

```python
def single_variable_groth(coeff_dict, varnum, beta)
```

Multiply ``sum_u coeff_u G_u^(beta)`` by the single variable ``x_varnum``
(Grothendieck Chevalley formula), via ``chevalley_x_k``.

<a id="schubmult.mult.groth.mult_poly_groth"></a>

#### mult\_poly\_groth

```python
def mult_poly_groth(coeff_dict, poly, var_x, beta)
```

Multiply ``sum_u coeff_u G_u^(beta)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable_groth``.

<a id="schubmult.mult.groth.groth_elem_sym_coeff"></a>

#### groth\_elem\_sym\_coeff

```python
def groth_elem_sym_coeff(k, u1, u2, vdiff, beta, length=None)
```

Coefficient of ``G_{u2}`` in ``E_{k - vdiff, k}(x; 0) G_{u1} = e_{k - vdiff}(x_1..x_k) G_{u1}``.

The ``y = z = 0`` specialization of ``groth_double._groth_elem_sym_frac``.  Sort the
window positions ``j <= k`` by the fate of ``u1(j)`` in ``u2``: *fixed*, *left*
(reappears at an earlier window position) or *out* (leaves the window or moves
right); with ``F = `fixed```, ``L = `left```, ``m = L + `out``` movers and
``d = l(u2) - l(u1)`` the closed form

    beta^(d - m) (-beta)^L E_{n - q, n}( (-)y_fixed, (-1/beta)^L ; z ),   n = F + L, q = vdiff,

collapses at ``y = z = 0`` (the fixed alphabet entries become ``0``) to

    beta^(d - m) (-beta)^(L - p) binom(L, p),   p = n - q,   0 <= p <= L,

and to ``0`` otherwise.  At ``beta = 0`` only ``d = m``, ``p = L`` survive, i.e. the
classical Pieri rule ``e_p(x_1..x_k) S_u = sum_{u ->_k w} S_w`` with ``p = k - d``
(there ``L = 0`` for chains with distinct lower indices, so ``p = F = k - d``).

``length`` overrides ``d``; the quantum kernel passes the quantum Bruhat chain length.

<a id="schubmult.mult.groth.grothmult_py"></a>

#### grothmult\_py

```python
def grothmult_py(perm_dict, v, beta=None)
```

Multiply (single) Grothendieck polynomials, mirroring ``schubmult_py``.

Computes the expansion of ``sum_u coeff_u G_u(x) * G_v(x)`` in the basis ``{G_w(x)}``
and returns it as ``{w: coeff_w}`` with coefficients polynomial in ``beta``.

Same method as ``groth_double.grothmult_double`` specialized to ``y = z = 0``: expand
``G_v`` into Schubert polynomials (``groth_elem_as_schub_dict``, via WC graphs) and
push each ``S_{v'}`` through the v-path kernel ``_groth_schub_vpath_mul``, whose
per-layer coefficients are the binomial closed form ``groth_elem_sym_coeff``.

