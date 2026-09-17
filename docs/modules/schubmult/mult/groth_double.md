<a id="schubmult.mult.groth_double"></a>

# schubmult.mult.groth\_double

K-theoretic Monk formula for double Grothendieck polynomials.

Implements Lenart--Postnikov, *Affine Weyl groups in K-theory and representation
theory* (arXiv:math/0309207), Corollary 8.2 (the :math:`K_T`-Monk formula), in
type :math:`A` and transported to the ``schubmult`` conventions.

Reconciling the conventions
---------------------------
``schubmult`` has no exponentials: it uses the multiplicative formal group law
``a (+) b = a + b + beta*a*b`` with formal inverse ``(-)b = -b/(1 + beta*b)``.
The double Grothendieck polynomial of a simple reflection is

.. math::

    \mathfrak{G}_{s_k}(x, y)
        = \frac{1}{\beta}\Bigl(\prod_{i=1}^{k}(1 + \beta x_i)(1 + \beta y_i) - 1\Bigr),

which at ``beta = -1`` is precisely the class ``1 - x^{w_0(omega_k)} e^{-omega_k}``
of Lemma 8.1(a).  The dictionary between the paper and this package is

* torus characters: ``x^{eps_i}  <->  1 + beta*y_i`` (so that ``(x^{eps_a - eps_b} - 1)/beta``
  is the root ``y_a (-) y_b``, matching ``DoubleGrothendieckRing.exp_root``);
* basis elements: ``[O_{X_{w_0 w}}]  <->  (-beta)^{l(w)} G_w``, since the paper
  indexes structure sheaves by the *dimension* of the Schubert variety whereas
  ``G_w`` has lowest term ``S_w`` of degree ``l(w)`` (codimension indexing).

Under ``u -> w_0 u`` the saturated *decreasing* chains of Corollary 8.2 become
saturated *increasing* chains, the reflections ``t_{ij}`` are unchanged, and the
character prefactor ``x^{nu(J)} = x^{w_0(omega_k) - u(omega_k)}`` (constant in
``J`` because ``omega_k`` is minuscule in type ``A``) becomes
``prod_{i<=k} (1 + beta*y_i)/(1 + beta*y_{u(i)})``.  Rescaling by
``(-beta)^{l(.)}`` turns the signs ``(-1)^{|J|}`` into powers ``beta^{|J|}`` and
yields

.. math::

    \mathfrak{G}_u(x, y)\,\mathfrak{G}_{s_k}(x, z)
        = \frac{1}{\beta}\Bigl(\Theta_u \sum_J \beta^{|J|}\,
          \mathfrak{G}_{u\,r_J}(x, y) - \mathfrak{G}_u(x, y)\Bigr),
    \qquad
    \Theta_u = \prod_{i=1}^{k}\frac{1 + \beta z_i}{1 + \beta y_{u(i)}},

the sum being over the subsets ``J`` of a reduced ``(-omega_k)``-chain of
reflections whose reflections build a saturated increasing Bruhat chain from
``u`` (the empty subset included).  The two secondary alphabets are handled by
``G_{s_k}(x, z) = C G_{s_k}(x, y) + (C - 1)/beta`` with
``C = prod_{i<=k}(1 + beta*z_i)/(1 + beta*y_i)``, which is exactly what turns
the ``y_i`` of ``x^{nu}`` into the ``z_i`` of ``Theta_u``.

By Corollary 15.4 of the same paper a reduced ``(-omega_k)``-chain of
reflections in type ``A_{n-1}`` is

    ``t_{1,n}, t_{1,n-1}, ..., t_{1,k+1}, t_{2,n}, ..., t_{2,k+1}, ..., t_{k,k+1}``.

<a id="schubmult.mult.groth_double.monk_chain"></a>

#### monk\_chain

```python
def monk_chain(k)
```

Reduced ``(-omega_k)``-chain of reflections in ``A_{n-1}``, as ``(i, j)``, ``i <= k < j``.

``omega_k = eps_1 + ... + eps_k``, so this is ``epsilon_chain`` on ``{1, ..., k}``:
every root ``alpha_{ij}`` with ``i, j <= k`` pairs to zero and drops out, and the
survivors all sit at level ``1``.  The order is ``i`` decreasing, then ``j`` decreasing.

<a id="schubmult.mult.groth_double.epsilon_chain"></a>

#### epsilon\_chain

```python
def epsilon_chain(positions, inverse=False, ambient_rank=None)
```

Reduced ``(-eps_A)``-chain of reflections in ``A_{n-1}``, ``A = positions``.

``positions`` is a single index or an iterable of them (repeats allowed), and
``eps_A = sum_{i in A} eps_i``.  With ``inverse=True`` the chain is for ``+eps_A``
instead, which is the weight of the inverse class ``prod_{i in A}(1 + beta*x_i)^{-1}``.

``(-omega_k)``-chains only see the roots ``alpha_{ij}`` with ``i <= k < j``, which is
why ``monk_chain`` multiplies by the whole product ``prod_{i<=k}(1 + beta*x_i)``.
Selecting an arbitrary set of variables needs ``eps_A`` instead, whose chain also
involves the roots ``alpha_{ik}`` with ``i < k``.

Built by Prop. 6.7: the reflections ``s_{alpha, m}`` separating the fundamental
alcove from ``A_{eps_A}``, ordered by the lexicographic key
``(lambda, alpha^vee)^{-1} (-m, (omega_1, alpha^vee), ..., (omega_{n-1}, alpha^vee))``.
Entries are ``(a, b, m)`` for the positive root ``alpha_{ab} = eps_a - eps_b``, ``a < b``;
``m > 0`` means ``b(r) = -alpha`` is negative and the step carries a sign in Thm 6.1.

Concatenating the individual ``(-eps_i)``-chains would also be legal (Prop. 12.2) but
only after translating the blocks, which shifts their levels; going through Prop. 6.7
avoids that and is reduced.  Note a single ``k`` gives ``(i, k, 0)`` for ``i < k`` and
``(k, j, 1)`` for ``j > k``, while for ``A = {1, ..., k}`` every ``alpha_{ij}`` with
``i, j <= k`` pairs to zero and drops out, leaving exactly ``monk_chain(k)``.
Flipping to ``inverse=True`` exchanges those two families.

<a id="schubmult.mult.groth_double.one_plus_beta_x_groth"></a>

#### one\_plus\_beta\_x\_groth

```python
def one_plus_beta_x_groth(coeff_dict,
                          positions,
                          var2=None,
                          beta=None,
                          inverse=False)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by ``prod_{i in positions} (1 + beta*x_i)``.

The one-pass Pieri rule of Theorem 6.1 at ``lambda = -eps_A``; see
``_one_plus_beta_x_terms`` for the coefficient.  ``inverse=True`` gives the inverse
operator ``prod_{i in positions} (1 + beta*x_i)^{-1}``, i.e. ``lambda = +eps_A``.

<a id="schubmult.mult.groth_double.single_variable_groth"></a>

#### single\_variable\_groth

```python
def single_variable_groth(coeff_dict, varnum, var2=None, beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by the single variable ``x_varnum``.

Returns ``{w: coeff_w}``.  This is ``_one_plus_beta_x_terms`` with the
identity subtracted off and ``beta`` divided out; the diagonal coefficient
collapses to ``(-) var2[u(varnum)] = -y/(1 + beta*y)``, the formal inverse of
``var2[u(varnum)]``, which is the localization of ``x_varnum`` at ``u``.

<a id="schubmult.mult.groth_double.elem_sym_perms_groth"></a>

#### elem\_sym\_perms\_groth

```python
def elem_sym_perms_groth(u, k)
```

K-theoretic analogue of ``elem_sym_perms``: ``{w: {d: multiplicity}}``.

Same recursion as ``elem_sym_perms(u, p, k)`` -- a step is any Bruhat cover
``w -> w t_{ij}`` with ``i <= k < j``, and ``j`` is required to weakly decrease along
the chain -- but with the ``p`` cut-off dropped, so chains of every length are
produced and the degree cut-off is left to the coefficient.

This is deliberately *not* a subset-of-a-fixed-chain enumeration.  A
``lambda``-chain imposes a total order on the transpositions, which loses covers such
as ``id < [1,3,2] < [2,3,1]`` (that needs ``t_{23}`` before ``t_{13}``); covers come
from arbitrary upward transpositions, and only ``j`` is constrained.

``d = l(w) - l(u)`` is the chain length.  A position ``i <= k`` may be stepped on more
than once, and distinct chains can land on the same ``w`` at the same ``d``, which is
the source of the K-theoretic multiplicities.

<a id="schubmult.mult.groth_double.grothmult_double_block"></a>

#### grothmult\_double\_block

```python
def grothmult_double_block(coeff_dict,
                           positions,
                           zvar=None,
                           var2=None,
                           beta=None,
                           fgl=True)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by a linear block over ``positions``:

    fgl=True  ->  prod_{i in A} (x_i (+) zvar),   x (+) z = x*(1 + beta*z) + z
    fgl=False ->  prod_{i in A} (x_i - zvar)

``positions`` is an arbitrary index set (repeats allowed), matching the ``index_list``
that ``pull_out_var`` produces, so this is the ``G``-basis analogue of the top-degree
mixed-variable block driving ``schubmult_double_alt`` / ``DoubleSchubertRing.elem_mul``.

``fgl=False`` is the plain ``beta = 0`` block ``(x_1 - z)(x_2 - z)...`` -- still a
perfectly good operator on the ``G`` basis, and the two are interchangeable via
``x (+) z = (1 + beta*z) * (x - (-)z)``, so either can be recovered from the other by
rescaling ``zvar``.

Computed by folding ``single_variable_groth`` one position at a time, using
``(a x_i + b) F = a (x_i F) + b F``.  That keeps every intermediate coefficient
polynomial in ``beta``; the one-pass alternative would expand
``prod_i ((1 + beta*x_i)(1 + beta*z) - 1) / beta**|A|`` by inclusion-exclusion over the
subsets of ``A`` (each term a ``one_plus_beta_x_groth`` call) and only cancel the
``beta^{-|A|}`` at the very end.

<a id="schubmult.mult.groth_double.groth_elem_sym_poly"></a>

#### groth\_elem\_sym\_poly

```python
def groth_elem_sym_poly(p, k, zvar, var_x, beta)
```

``E_p^beta(x_1..x_k; z) = e_p(x_1 (+) z, ..., x_k (+) z)``, ``x (+) z = x(1 + beta*z) + z``.

The double Grothendieck elementary symmetric: ``p == k`` gives
``(x_1(1 + beta*z) + z) ... (x_k(1 + beta*z) + z)`` and ``beta == 0`` gives the
factorial elementary symmetric ``elem_sym_poly(p, k, x, [-z])``.

<a id="schubmult.mult.groth_double.grothmult_double_pieri"></a>

#### grothmult\_double\_pieri

```python
def grothmult_double_pieri(coeff_dict,
                           p,
                           k,
                           zvar=None,
                           var_x=None,
                           var2=None,
                           beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by ``groth_elem_sym_poly(p, k, zvar, var_x, beta)``.

Exact, by folding the ``e_p`` DP over coefficient dicts one shifted variable
``x_i (+) zvar`` at a time (``single_variable_groth`` per step) -- no symbolic
expansion.  A closed-form Pieri rule in the style of ``dom_groth`` -- paths from
``elem_sym_perms_groth`` plus an ``elem_sym_poly`` in the localizations -- is *not*
implemented: grading the paths by ``beta^{d - m}`` with ``m`` the number of moved
positions and taking ``elem_sym_poly`` over the untouched ones is wrong already at
``p == k``.  The non-equivariant rule ``groth_pieri_mul`` grades instead by
``beta^{d - (number of marked steps)}`` with the multiplicity counting admissible
markings of the chain (``elem_sym_chains_groth``), so the equivariant coefficient
presumably needs that marking data rather than the moved/untouched split.

``var_x`` is unused (kept for signature compatibility).

<a id="schubmult.mult.groth_double.grothmult_double_top"></a>

#### grothmult\_double\_top

```python
def grothmult_double_top(coeff_dict, k, zvar=None, var2=None, beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by the top linear block ``prod_{i=1}^{k}(x_i + zvar)``.

Closed positive Molev--Sagan Pieri rule (conjectural; verified exhaustively on
``S_4`` and sampled through ``S_6``, ``k <= 5``, against
``grothmult_double_block(..., zvar=-zvar, fgl=False)``):

.. math::

    \prod_{i=1}^{k}(x_i + z)\,\mathfrak{G}_u(x; y)
        = \sum_{w} \beta^{\,d - k + |Q|} \Bigl(\prod_{i=1}^{k} f_i\Bigr)\,
          \mathfrak{G}_w(x; y),
    \qquad d = \ell(w) - \ell(u),

where the factor ``f_i`` depends on the fate of the window value ``u(i)``:

* ``u(i) = w(i)`` (the set ``Q``):  ``(z(1 + beta*y_{u(i)}) - y_{u(i)}) / (1 + beta*y_{u(i)})``,
  i.e. ``z (+) (-)y_{u(i)}``, the K-theoretic analogue of ``z - y_{u(i)}``;
* ``u(i)`` stays in the window but moves left:  ``1 - beta*zvar``;
* ``u(i)`` exits the window or moves right within it:  ``1/(1 + beta*y_{u(i)})``.

The sum runs over the marked-chain K-Pieri support (``_top_block_support``).
At ``beta = 0`` this collapses to the ``p = k`` Pieri formula for double
Schubert polynomials [Samuel, Theorem 7.1]:
``S_u(x;y) prod(x_i - z) = sum_{u ->_k w} prod_{i in Q}(y_{u(i)} - z) S_w(x;y)``
with ``z -> -z``.

<a id="schubmult.mult.groth_double.mult_poly_groth_double"></a>

#### mult\_poly\_groth\_double

```python
def mult_poly_groth_double(coeff_dict,
                           poly,
                           var_x=None,
                           var_y=None,
                           beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var_y)`` by an arbitrary polynomial in ``var_x``.

Mirrors ``mult_poly_double``; the leaves of the ``Add``/``Mul``/``Pow`` recursion
are handled by ``single_variable_groth``.

<a id="schubmult.mult.groth_double.dgroth_to_dschub"></a>

#### dgroth\_to\_dschub

```python
def dgroth_to_dschub(v, var3, beta=None)
```

Expand ``G_v(x, var3)`` in double Schubert polynomials: ``{v': coeff}``.

``sum_{v'} coeff_{v'} S_{v'}(x, var3) = G_v(x, var3)`` with coefficients in
``var3`` and ``beta``.  Exact but slow; delegates to ``grothendieck_poly``
with ``keep_as_schub=True``.

<a id="schubmult.mult.groth_double.groth_elem_sym_func"></a>

#### groth\_elem\_sym\_func

```python
def groth_elem_sym_func(k, i, u1, u2, v1, v2, vdiff, varl1, varl2, beta)
```

Expression form of ``_groth_elem_sym_frac``; see there for the rule.

<a id="schubmult.mult.groth_double.grothmult_double"></a>

#### grothmult\_double

```python
def grothmult_double(perm_dict, v, var2=None, var3=None, beta=None)
```

Multiply double Grothendieck polynomials, mirroring ``schubmult_double``.

Computes the expansion of ``sum_u coeff_u G_u(x, var2) * G_v(x, var3)`` in
the basis ``{G_w(x, var2)}`` and returns it as ``{w: coeff_w}``.

``v = s_k`` uses the verified chain formula of Corollary 8.2, and
``max_descent == 1`` folds that column by column.  General ``v`` goes through
``dgroth_to_dschub`` (exact, slow) and the conjectural vpath kernel
``_groth_schub_vpath_mul``, one run per double Schubert ``S_{v'}`` in the
expansion of ``G_v``.

The chain rank is inferred from the current permutation and selected positions.

