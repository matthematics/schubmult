<a id="schubmult.mult.groth_quantum_double"></a>

# schubmult.mult.groth\_quantum\_double

Quantum double Grothendieck multiplication: the Molev--Sagan machinery for ``QK_T(Fl_n)``.

Conventions
-----------
The quantum double Grothendieck polynomials are the Lenart--Maeno quantizations of the double
Grothendieck polynomials in the ``x`` alphabet, with ``beta`` and the secondary alphabet treated
as scalars.  Write ``X_i = 1 + beta*x_i`` (the K-theoretic line bundle variables; ``1 - x_i`` at
``beta = -1``).  The quantization map ``Q`` is the linear map that is multiplicative on the
standard elementary monomials ``prod_j e_{i_j}(X_1..X_j)``, ``i_j <= j``, and sends

    e_l(X_1..X_k)  ->  F^k_l(X) = sum_{J in [k], |J| = l} prod_{j in J, j+1 not in J} (1 - Q_j) prod_{j in J} X_j,
    Q_j = beta^2 q_j

(``lm_quantize``).  ``G^q_v(x; y) := Q(G_v(x; y))``; at ``beta = -1``, ``Q_j`` these are the
Lenart--Maeno polynomials that represent the Schubert classes of ``QK_T(Fl_n)`` in the
Maeno--Naito--Sagaki presentation (arXiv:2302.09485, 2305.17685), with ``e^{-eps_j} = 1 - y_j``.
The normalization ``Q_j = beta^2 q_j`` makes ``deg Q_j = 0`` (``deg beta = -1``, ``deg q_j = 2``)
and gives ``G^q_v(x; y)|_{beta = 0} = S^q_v(x; -y)`` (quantum double Schubert).  Note this is not
the Fomin--Gelfand--Postnikov quantization of the ``x`` alphabet: already
``G^q_{21}(x; y) = x_1 (+) y_1 - beta q_1 (1 + beta x_1)(1 + beta y_1)``.

Every function here returns ``{w: coeff}`` meaning ``sum_w coeff_w G^q_w(x; var2)``.  The products
are polynomial identities in ``Z[beta, q, y, z][x]``: the ``G^q_w`` are stable in ``n`` and span
the same filtered pieces as the classical ``G_w``, so no quotient by the quantum ideal is needed.

The rule
--------
Everything is read off the equivariant quantum top block

    Q( prod_{i<=k} (x_i + z) ) G^q_u(x; y)
        = sum_w beta^{len - m} q^D prod_{Fix}(z + (-)y_a) (1 - beta z)^{|Left|}
              prod_{Out} (1 + beta y_a)^{-1}  G^q_w(x; y),

where ``w`` runs over the endpoints of the Naito--Sagaki ``k``-Pieri chains
(arXiv:2211.01578) in the quantum Bruhat graph starting at ``u``, ``len`` is the length of the
chain and ``q^D`` its quantum weight, and (Fix, Left, Out) sorts the window values ``u(i)``,
``i <= k``, exactly as in ``grothmult_double_top`` (``m = |Left| + |Out|``).  Equivalently the
prefactor is ``beta^{l(w) - l(u) - m} Q^D``.  The pair ``(len, D)`` is an invariant of
``(u, w, k)`` in every case computed so far; the code raises if that ever fails.  The rule is
conjectural: it is verified in ``QK_T(Fl_3)`` symbolically and in ``QK_T(Fl_4)`` under random
specializations against the Maeno--Naito--Sagaki presentation
(``_lscripts/qk_equivariant_oracle.py``), and as a polynomial identity by
``_lscripts/qgroth_kernel_check.py``.  Its specializations are theorems or existing kernels:
``q = 0`` is ``grothmult_double_top``, ``beta = 0`` is the ``p = k`` quantum Pieri rule of
``schubmult_q_double``, and ``y = 0``, ``beta = -1`` is Naito--Sagaki's quantum K Pieri theorem.

Molev--Sagan assembly
---------------------
``S_{v'}(x; z)`` is a sum over strict-theta v-paths of products of factorial elementary symmetric
polynomials ``E_{p,k}(x; z)`` with strictly decreasing ``k``.  Each ``E_{p,k}`` is symmetric in
``x_1..x_k`` and of degree ``<= 1`` in each variable, hence a combination of the
``e_l(X_1..X_k)``, so the product is a combination of standard elementary monomials and
``Q`` is multiplicative on it: ``Q(S_{v'}(x; z))`` is the same v-path sum with ``Q(E_{p,k})``.
Since ``E_{k-q,k}(x; z_1..z_{q+1}) = (-1)^q d^z_q ... d^z_1 prod_{i<=k}(x_i - z_1)`` and ``Q``
commutes with the ``z``-divided differences, the coefficient of ``G^q_{u2}`` in
``Q(E_{k-q,k}(x; z)) G^q_{u1}(x; y)`` is ``_groth_elem_sym_frac`` with ``l(u2) - l(u1)`` replaced
by the chain length and the quantum weight ``q^D`` attached.  Chaining the layers of the v-path
recursion exactly as ``schubmult_q_double`` does gives ``G^q_u(x; y) Q(S_{v'}(x; z))``, and
``dgroth_to_dschub`` (``G_v = sum c_{v'} S_{v'}``, ``Q`` linear over ``z, beta``) turns that into
``G^q_u(x; y) G^q_v(x; z)``.

<a id="schubmult.mult.groth_quantum_double.quantum_pieri_chains"></a>

#### quantum\_pieri\_chains

```python
@cache
def quantum_pieri_chains(u, k)
```

Endpoints of the Naito--Sagaki ``k``-Pieri chains from ``u`` in the quantum Bruhat graph.

A ``k``-Pieri chain is a path ``u = w_0 -> w_1 -> ... -> w_r`` with edges
``w -> w t_{ab}``, ``a <= k < b``, that are either Bruhat covers or quantum edges
(``l`` drops by ``2(b - a) - 1``, weight ``q_a ... q_{b-1}``), whose labels ``(a, b)`` are
distinct, have ``b`` weakly decreasing, and satisfy: whenever a label repeats an earlier
lower index ``a``, the next label is larger in ``_ns_prec``.  Every such chain admits a
Naito--Sagaki marking, so these endpoints are exactly the support of the quantum top block
(at ``q = 0`` they reduce to ``_top_block_support``).

Returns ``{w: (length, D)}`` with ``D`` the tuple of ``q`` exponents (``D[j - 1]`` is the
exponent of ``q_j``).  The empty chain contributes ``u: (0, 0)``.  Raises ``ValueError``
if two chains to the same ``w`` disagree on ``(length, D)``, which would leave the rule
undefined.

<a id="schubmult.mult.groth_quantum_double.grothmult_q_double_top"></a>

#### grothmult\_q\_double\_top

```python
def grothmult_q_double_top(coeff_dict,
                           k,
                           zvar=None,
                           var2=None,
                           beta=None,
                           q_var=None)
```

Multiply ``sum_u coeff_u G^q_u(x, var2)`` by the quantized top block ``Q(prod_{i<=k}(x_i + zvar))``.

The multiplier is ``groth_elem_sym_poly_q(k, k, zvar, x, beta, q_var, fgl=False)``, i.e.
``beta^{-k} sum_l (-(1 - beta z))^{k-l} F^k_l(X)``.  The coefficient of ``G^q_w`` is
``_top_block_coeff`` evaluated with the quantum chain length, times the quantum weight
``q^D``; see the module docstring.

<a id="schubmult.mult.groth_quantum_double.groth_elem_sym_poly_q"></a>

#### groth\_elem\_sym\_poly\_q

```python
def groth_elem_sym_poly_q(p, k, zvar, var_x, beta, q_var=None, fgl=True)
```

Quantization of ``groth_elem_sym_poly``: ``Q(e_p(x_1 (+) z, ..., x_k (+) z))``.

``fgl=False`` quantizes the plain ``e_p(x_1 + z, ..., x_k + z)`` instead, whose ``p = k``
case is the top block of ``grothmult_q_double_top``.  Both are symmetric in ``x_1..x_k`` of
degree ``<= 1`` in each variable, so ``lm_quantize`` with ``k + 1`` slots applies.

<a id="schubmult.mult.groth_quantum_double.grothmult_q_double_pieri"></a>

#### grothmult\_q\_double\_pieri

```python
def grothmult_q_double_pieri(coeff_dict,
                             p,
                             k,
                             zvar=None,
                             var2=None,
                             beta=None,
                             q_var=None,
                             fgl=True)
```

Multiply ``sum_u coeff_u G^q_u(x, var2)`` by ``groth_elem_sym_poly_q(p, k, zvar, x, beta, q_var, fgl)``.

Closed form from the top block (see ``_q_pieri_coeff``); ``fgl=False`` multiplies by the
quantization of ``e_p(x_1 + z, ..., x_k + z)``.  At ``q = 0`` this agrees with the exact
fold ``grothmult_double_pieri``.

<a id="schubmult.mult.groth_quantum_double.grothmult_q_double"></a>

#### grothmult\_q\_double

```python
def grothmult_q_double(perm_dict,
                       v,
                       var2=None,
                       var3=None,
                       beta=None,
                       q_var=None)
```

Multiply quantum double Grothendieck polynomials, mirroring ``schubmult_q_double``.

Returns the expansion of ``sum_u coeff_u G^q_u(x, var2) * G^q_v(x, var3)`` in the basis
``{G^q_w(x, var2)}`` as ``{w: coeff_w}``, coefficients rational in ``var2`` (denominators
are products of ``1 + beta*var2[a]``) and polynomial in ``var3``, ``beta``, ``q``.

``G^q_v(x, var3)`` is expanded through ``dgroth_to_dschub`` and one run of the quantum
v-path kernel ``_qgroth_schub_vpath_mul`` per double Schubert term.  At ``q = 0`` this is
``grothmult_double``; at ``beta = 0`` it is ``schubmult_q_double``.

<a id="schubmult.mult.groth_quantum_double.grothmult_q_double_dict"></a>

#### grothmult\_q\_double\_dict

```python
def grothmult_q_double_dict(perm_dict1,
                            perm_dict2,
                            var2=None,
                            var3=None,
                            beta=None,
                            q_var=None)
```

Product of two coefficient dicts: ``sum_v coeff2_v grothmult_q_double(perm_dict1, v, ...)``.

<a id="schubmult.mult.groth_quantum_double.quantum_elem_sym"></a>

#### quantum\_elem\_sym

```python
def quantum_elem_sym(l, k, var_x, beta, q_var=None)
```

``F^k_l(X) = sum_{J in [k], |J| = l} prod_{j in J, j+1 not in J} (1 - beta^2 q_j) prod_{j in J} X_j``, ``X_j = 1 + beta*x_j``.

The Lenart--Maeno quantization of ``e_l(X_1..X_k)``; ``beta = -1`` gives the ``F^k_l`` of
Maeno--Naito--Sagaki with ``Q_j = q_j``.  No ``1 - Q_N := 1`` convention is applied: that
belongs to the defining ideal of ``QK_T(Fl_N)``, not to the polynomials.

<a id="schubmult.mult.groth_quantum_double.lm_quantize"></a>

#### lm\_quantize

```python
def lm_quantize(poly, N, var_x, beta, q_var=None)
```

Lenart--Maeno quantization of a polynomial in ``x_1..x_{N-1}`` of degree ``<= N - i`` in ``x_i``.

Expands ``poly`` in the standard elementary monomials ``prod_j e_{i_j}(X_1..X_j)`` of
``X_i = 1 + beta*x_i`` (a basis of that span; ``_sem_basis``) and replaces each
``e_{i_j}(X_1..X_j)`` by ``quantum_elem_sym(i_j, j)``.  Linear over everything but ``x``,
stable in ``N`` (``poly`` may use fewer variables), and ``Q_j = beta^2 q_j`` so the result
is polynomial in ``beta``.  Raises if ``poly`` is not in the span.

<a id="schubmult.mult.groth_quantum_double.qgroth_poly"></a>

#### qgroth\_poly

```python
def qgroth_poly(v, var_x=None, var_y=None, beta=None, q_var=None)
```

The quantum double Grothendieck polynomial ``G^q_v(var_x; var_y)`` as an explicit expression.

``lm_quantize`` applied to ``grothendieck_poly(v)`` with ``N = len(v)`` slots.  At
``beta = -1`` this is the Maeno--Naito--Sagaki ``G^Q_v(x, y)`` with ``Q_j = q_j``.  Slow
(symbolic); intended for verification.

