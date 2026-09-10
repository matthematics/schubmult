r"""K-theoretic Monk formula for double Grothendieck polynomials.

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
"""

from fractions import Fraction

from schubmult.abc import beta as _default_beta
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import Add, Mul, Pow, S, prod, sympify
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet_base
from schubmult.utils.perm_utils import add_perm_dict

__all__ = [
    "dgroth_to_dschub",
    "elem_sym_perms_groth",
    "epsilon_chain",
    "groth_elem_sym_func",
    "groth_elem_sym_poly",
    "grothmult_double",
    "grothmult_double_block",
    "grothmult_double_pieri",
    "grothmult_double_top",
    "monk_chain",
    "mult_poly_groth_double",
    "one_plus_beta_x_groth",
    "single_variable_groth",
]


def monk_chain(k):
    """Reduced ``(-omega_k)``-chain of reflections in ``A_{n-1}``, as ``(i, j)``, ``i <= k < j``.

    ``omega_k = eps_1 + ... + eps_k``, so this is ``epsilon_chain`` on ``{1, ..., k}``:
    every root ``alpha_{ij}`` with ``i, j <= k`` pairs to zero and drops out, and the
    survivors all sit at level ``1``.  The order is ``i`` decreasing, then ``j`` decreasing.
    """
    return tuple((a, b) for a, b, _ in epsilon_chain(tuple(range(1, k + 1))))


def epsilon_chain(positions, inverse=False, ambient_rank=None):
    r"""Reduced ``(-eps_A)``-chain of reflections in ``A_{n-1}``, ``A = positions``.

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
    """
    if isinstance(positions, int):
        positions = (positions,)
    step = 1 if inverse else -1
    if ambient_rank is None:
        ambient_rank = max(positions) + len(positions)
    weight = [0] * ambient_rank
    for k in positions:
        weight[k - 1] += step

    entries = []
    for a in range(1, ambient_rank + 1):
        for b in range(a + 1, ambient_rank + 1):
            pairing = weight[a - 1] - weight[b - 1]
            if pairing > 0:
                levels = range(0, -pairing, -1)
            elif pairing < 0:
                levels = range(1, -pairing + 1)
            else:
                continue
            omegas = tuple((1 if a <= r else 0) - (1 if b <= r else 0) for r in range(1, ambient_rank))
            for m in levels:
                key = tuple(Fraction(value, pairing) for value in (-m, *omegas))
                entries.append((key, a, b, m))
    entries.sort(key=lambda entry: entry[0])
    return tuple((a, b, m) for _, a, b, m in entries)


def _one_plus_beta_x_terms(u, positions, var2, beta, inverse=False, beta_offset=0):
    r"""Expand ``prod_{i in A} (1 + beta*x_i)^{-1 if inverse else 1} * G_u(x, var2)``.

    ``prod_{i in A}(1 + beta*x_i)`` is the class ``e^{-eps_A}``, so this is Theorem 6.1
    at ``lambda = -eps_A``, in one pass over ``epsilon_chain(positions)``.  Writing
    ``J`` for a subset of that chain whose reflections form a saturated increasing
    Bruhat chain ``u = w_0 < w_1 < ... < w_s = w``, the transported coefficient is

        (-1)^{n(J)} (-beta)^{|J|} * prod_{i in A} (1 + beta*var2[w(i)])^{-1}
            * prod_{steps} ( (1 + beta*var2[w_j(b)]) / (1 + beta*var2[w_j(a)]) )^{level},

    where ``w_j`` is the permutation just before the step and ``n(J)`` counts the steps
    with a negative root, i.e. those of positive level.  Distinct ``J`` can land on the
    same ``w``, which is where the K-theoretic multiplicities come from, and a position
    may be stepped on more than once -- both impossible in the ``beta = 0`` rule.

    ``inverse=True`` is the class ``e^{+eps_A}``: the chain is taken for ``+eps_A`` and,
    since ``-mu = w(-lambda) + u(tau)`` flips with ``lambda``, the terminal factor moves
    to the numerator.  The per-step translation factor is unchanged.
    """
    chain = epsilon_chain(positions, inverse=inverse, ambient_rank=_rank(u, positions))
    exponent = 1 if inverse else -1
    terms = {}

    def terminal(w):
        value = S.One
        for i in positions:
            value *= S.One + beta * var2[w[i - 1]]
        return value**exponent

    def walk(start, w, character, sign, depth):
        # beta_offset strips beta factors structurally (chains of depth > 0 all
        # carry beta**depth); the depth-0 diagonal term is skipped when offset
        if depth or not beta_offset:
            terms[w] = terms.get(w, S.Zero) + sign * beta ** (depth + beta_offset) * character * terminal(w)
        for index in range(start, len(chain)):
            a, b, level = chain[index]
            stepped = w.swap(a - 1, b - 1)
            if stepped.inv != w.inv + 1:
                continue
            # s_{alpha_ab, level} translates the weight by level * w(alpha_ab).
            factor = ((S.One + beta * var2[w[b - 1]]) / (S.One + beta * var2[w[a - 1]])) ** level
            walk(index + 1, stepped, character * factor, sign if level > 0 else -sign, depth + 1)

    walk(0, u, S.One, S.One, 0)
    return terms


def _rank(u, positions):
    return max(len(u), max(positions) + 1) + len(positions)


def one_plus_beta_x_groth(coeff_dict, positions, var2=None, beta=None, inverse=False):
    r"""Multiply ``sum_u coeff_u G_u(x, var2)`` by ``prod_{i in positions} (1 + beta*x_i)``.

    The one-pass Pieri rule of Theorem 6.1 at ``lambda = -eps_A``; see
    ``_one_plus_beta_x_terms`` for the coefficient.  ``inverse=True`` gives the inverse
    operator ``prod_{i in positions} (1 + beta*x_i)^{-1}``, i.e. ``lambda = +eps_A``.
    """
    if beta is None:
        beta = _default_beta
    var2 = _genset(var2)
    positions = (positions,) if isinstance(positions, int) else tuple(positions)

    ret = {}
    for u, val in coeff_dict.items():
        u = Permutation(u)
        for w, coeff in _one_plus_beta_x_terms(u, positions, var2, beta, inverse=inverse).items():
            ret[w] = ret.get(w, S.Zero) + val * coeff
    return ret


def single_variable_groth(coeff_dict, varnum, var2=None, beta=None):
    r"""Multiply ``sum_u coeff_u G_u(x, var2)`` by the single variable ``x_varnum``.

    Returns ``{w: coeff_w}``.  This is ``_one_plus_beta_x_terms`` with the
    identity subtracted off and ``beta`` divided out; the diagonal coefficient
    collapses to ``(-) var2[u(varnum)] = -y/(1 + beta*y)``, the formal inverse of
    ``var2[u(varnum)]``, which is the localization of ``x_varnum`` at ``u``.
    """
    if beta is None:
        beta = _default_beta
    var2 = _genset(var2)
    k = varnum
    positions = (k,)

    ret = {}
    for u, val in coeff_dict.items():
        u = Permutation(u)
        y = var2[u[k - 1]]
        ret[u] = ret.get(u, S.Zero) + val * (-y) / (S.One + beta * y)
        # beta_offset=-1 divides out the beta structurally, no cancel() needed
        for w, coeff in _one_plus_beta_x_terms(u, positions, var2, beta, beta_offset=-1).items():
            ret[w] = ret.get(w, S.Zero) + val * coeff
    return ret


def elem_sym_perms_groth(u, k):
    r"""K-theoretic analogue of ``elem_sym_perms``: ``{w: {d: multiplicity}}``.

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
    """
    out = {}

    def walk(w, last_b, d):
        if d:
            counts = out.setdefault(w, {})
            counts[d] = counts.get(d, 0) + 1
        for b in range(last_b, k, -1):
            for a in range(1, k + 1):
                stepped = w.swap(a - 1, b - 1)
                if stepped.inv == w.inv + 1:
                    walk(stepped, b, d + 1)

    walk(u, _rank(u, tuple(range(1, k + 1))), 0)
    return out


def grothmult_double_block(coeff_dict, positions, zvar=None, var2=None, beta=None, fgl=True):
    r"""Multiply ``sum_u coeff_u G_u(x, var2)`` by a linear block over ``positions``:

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
    """
    if beta is None:
        beta = _default_beta
    if zvar is None:
        zvar = S.Zero
    var2 = _genset(var2)
    positions = (positions,) if isinstance(positions, int) else tuple(positions)

    scale, shift = (S.One + beta * zvar, zvar) if fgl else (S.One, -zvar)

    ret = {Permutation(key): value for key, value in coeff_dict.items()}
    for i in positions:
        stepped = {w: scale * coeff for w, coeff in single_variable_groth(ret, i, var2, beta).items()}
        for w, coeff in ret.items():
            stepped[w] = stepped.get(w, S.Zero) + shift * coeff
        ret = stepped
    return ret


def groth_elem_sym_poly(p, k, zvar, var_x, beta):
    """``E_p^beta(x_1..x_k; z) = e_p(x_1 (+) z, ..., x_k (+) z)``, ``x (+) z = x(1 + beta*z) + z``.

    The double Grothendieck elementary symmetric: ``p == k`` gives
    ``(x_1(1 + beta*z) + z) ... (x_k(1 + beta*z) + z)`` and ``beta == 0`` gives the
    factorial elementary symmetric ``elem_sym_poly(p, k, x, [-z])``.
    """
    acc = [S.One] + [S.Zero] * p
    for i in range(1, k + 1):
        shifted = var_x[i] * (S.One + beta * zvar) + zvar
        for r in range(p, 0, -1):
            acc[r] = acc[r] + acc[r - 1] * shifted
    return acc[p]


def grothmult_double_pieri(coeff_dict, p, k, zvar=None, var_x=None, var2=None, beta=None):  # noqa: ARG001
    r"""Multiply ``sum_u coeff_u G_u(x, var2)`` by ``groth_elem_sym_poly(p, k, zvar, var_x, beta)``.

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
    """
    if beta is None:
        beta = _default_beta
    if zvar is None:
        zvar = S.Zero
    var2 = _genset(var2)

    # acc[r] = e_r(x_1 (+) z, ..., x_i (+) z) * input, folded as coefficient dicts
    acc = [{Permutation(key): value for key, value in coeff_dict.items()}] + [{} for _ in range(p)]
    scale = S.One + beta * zvar
    for i in range(1, k + 1):
        for r in range(min(p, i), 0, -1):
            if not acc[r - 1]:
                continue
            term = {w: scale * c for w, c in single_variable_groth(acc[r - 1], i, var2, beta).items()}
            if zvar != S.Zero:
                for w, c in acc[r - 1].items():
                    term[w] = term.get(w, S.Zero) + zvar * c
            acc[r] = add_perm_dict(acc[r], term)
    return acc[p]


def _top_block_support(u, k):
    """Support of ``prod_{i<=k}(x_i + z) G_u``: endpoints of the marked K-Pieri chains.

    Equals ``union_{p=0..k} supp(e_p^beta(x_1..x_k) G_u)``.  A marked chain from
    ``elem_sym_chains_groth`` contributes to degree ``p`` iff
    ``fff <= p <= fff + (len - 1 - fff - ppp)`` with ``fff``/``ppp`` the forced
    marks/unmarks, so the union over ``0 <= p <= k`` is nonempty iff ``fff <= k``.
    This is strictly smaller than the ``elem_sym_perms_groth`` endpoint set: a
    ``j``-weakly-decreasing cover chain need not admit a valid marking (e.g.
    ``u = [1,2,4,3]``, ``k = 2``, ``w = [3,4,1,2]``).
    """
    from schubmult.utils.schub_lib import elem_sym_chains_groth

    support = set()
    for perms, markings in elem_sym_chains_groth(u, 0, k):
        if sum(1 for m in markings if m == 1) <= k:
            support.add(perms[-1])
    return support


def _top_block_coeff(u, w, k, zvar, var2, beta):
    """Coefficient of ``G_w`` in ``prod_{i<=k}(x_i + zvar) G_u``: one factor per window position."""
    window = [w[i] for i in range(k)]
    fixed = 0
    value = S.One
    for i in range(k):
        v = u[i]
        y = var2[v]
        if window[i] == v:
            fixed += 1
            value *= (zvar * (S.One + beta * y) - y) / (S.One + beta * y)
        elif v in window and window.index(v) < i:
            value *= S.One - beta * zvar
        else:
            value *= S.One / (S.One + beta * y)
    power = (w.inv - u.inv) - (k - fixed)
    if power < 0:
        raise ValueError(f"negative beta power on support element: u={list(u)}, w={list(w)}, k={k}")
    return beta**power * value


def grothmult_double_top(coeff_dict, k, zvar=None, var2=None, beta=None):
    r"""Multiply ``sum_u coeff_u G_u(x, var2)`` by the top linear block ``prod_{i=1}^{k}(x_i + zvar)``.

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
    """
    if beta is None:
        beta = _default_beta
    if zvar is None:
        zvar = S.Zero
    var2 = _genset(var2)

    ret = {}
    for u, val in coeff_dict.items():
        u = Permutation(u)
        if k == 0:
            ret[u] = ret.get(u, S.Zero) + val
            continue
        for w in _top_block_support(u, k):
            ret[w] = ret.get(w, S.Zero) + val * _top_block_coeff(u, w, k, zvar, var2, beta)
    return ret


def mult_poly_groth_double(coeff_dict, poly, var_x=None, var_y=None, beta=None):
    """Multiply ``sum_u coeff_u G_u(x, var_y)`` by an arbitrary polynomial in ``var_x``.

    Mirrors ``mult_poly_double``; the leaves of the ``Add``/``Mul``/``Pow`` recursion
    are handled by ``single_variable_groth``.
    """
    var_x = _genset(var_x)
    var_y = _genset(var_y)
    index = var_x.index(poly)
    if index != -1:
        return single_variable_groth(coeff_dict, index, var_y, beta)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for arg in poly.args:
            ret = mult_poly_groth_double(ret, arg, var_x, var_y, beta)
        return ret
    if isinstance(poly, Pow):
        base, exponent = poly.args
        ret = coeff_dict
        for _ in range(int(exponent)):
            ret = mult_poly_groth_double(ret, base, var_x, var_y, beta)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for arg in poly.args:
            ret = add_perm_dict(ret, mult_poly_groth_double(coeff_dict, arg, var_x, var_y, beta))
        return ret
    return {perm: poly * coeff for perm, coeff in coeff_dict.items()}


def _chain_sums(u, k, beta):
    """``{w: sum_J beta**(|J| - 1)}`` over nonempty ``J`` with ``u r_J = w``.

    ``J`` runs over subsets of ``monk_chain(k)`` whose reflections, applied in
    chain order, form a saturated increasing chain in Bruhat order from ``u``.
    """
    chain = monk_chain(k)
    sums = {}

    def walk(start, w, size):
        if size:
            sums[w] = sums.get(w, S.Zero) + beta ** (size - 1)
        for index in range(start, len(chain)):
            a, b = chain[index]
            stepped = w.swap(a - 1, b - 1)
            if stepped.inv == w.inv + 1:
                walk(index + 1, stepped, size + 1)

    walk(0, u, 0)
    return sums


def _genset(var):
    if var is None or isinstance(var, GeneratingSet_base):
        return var
    return CustomGeneratingSet(var)


def dgroth_to_dschub(v, var3, beta=None):
    """Expand ``G_v(x, var3)`` in double Schubert polynomials: ``{v': coeff}``.

    ``sum_{v'} coeff_{v'} S_{v'}(x, var3) = G_v(x, var3)`` with coefficients in
    ``var3`` and ``beta``.  Exact but slow; delegates to ``grothendieck_poly``
    with ``keep_as_schub=True``.
    """
    from schubmult.symbolic.poly.schub_poly import grothendieck_poly
    from schubmult.symbolic.poly.variables import GeneratingSet

    if beta is None:
        beta = _default_beta
    var3 = _genset(var3)
    elem = grothendieck_poly(Permutation(v), GeneratingSet("x"), var3, beta, keep_as_schub=True)
    return {Permutation(key): value for key, value in elem.items() if value != S.Zero}


def _mul_linear(coeffs, c0, c1):
    """Multiply the polynomial ``sum_m coeffs[m] t^m`` by ``c0 + c1*t``."""
    out = [S.Zero] * (len(coeffs) + 1)
    for m, c in enumerate(coeffs):
        out[m] += c * c0
        out[m + 1] += c * c1
    return out


def _complete_homog(p, vrs):
    """Complete homogeneous symmetric polynomial ``h_p`` of ``vrs``; ``h_{<0} = 0``."""
    if p < 0:
        return S.Zero
    acc = [S.One] + [S.Zero] * p
    for v in vrs:
        for q in range(1, p + 1):
            acc[q] = acc[q] + v * acc[q - 1]
    return acc[p]


def _frac_mul(f1, f2):
    """Multiply two flat fractions ``(numer, {atom: exp})``; denominators multiply by adding exponents."""
    n1, d1 = f1
    n2, d2 = f2
    d = dict(d1)
    for a, e in d2.items():
        d[a] = d.get(a, 0) + e
    return (n1 * n2, d)


def _frac_add(f1, f2, varl1, beta):
    """Add two flat fractions over the LCM denominator.

    Each side is only scaled by the *missing* atoms ``(1 + beta*varl1[a])``, so the
    result stays a single-level fraction with a canonical atom-exponent denominator.
    The numerator is left unexpanded (small, but possibly a hidden zero).
    """
    if f1 is None:
        return f2
    n1, d1 = f1
    n2, d2 = f2
    d = {a: max(d1.get(a, 0), d2.get(a, 0)) for a in set(d1) | set(d2)}
    m1 = prod([(S.One + beta * varl1[a]) ** (d[a] - d1.get(a, 0)) for a in d])
    m2 = prod([(S.One + beta * varl1[a]) ** (d[a] - d2.get(a, 0)) for a in d])
    return (n1 * m1 + n2 * m2, d)


def _frac_to_expr(f, varl1, beta):
    n, d = f
    if n == S.Zero:
        return S.Zero
    return n / prod([(S.One + beta * varl1[a]) ** e for a, e in d.items() if e])


def groth_elem_sym_func(k, i, u1, u2, v1, v2, vdiff, varl1, varl2, beta):
    """Expression form of ``_groth_elem_sym_frac``; see there for the rule."""
    return _frac_to_expr(_groth_elem_sym_frac(k, i, u1, u2, v1, v2, vdiff, varl1, varl2, beta), varl1, beta)


def _groth_elem_sym_frac(k, i, u1, u2, v1, v2, vdiff, varl1, varl2, beta):
    r"""K-analogue of ``elem_sym_func`` for the vpath iteration.

    Coefficient of ``G_{u2}(x, varl1)`` contributed when layer ``i`` (block size
    ``k``) multiplies ``G_{u1}(x, varl1)`` while the v-path steps ``v1 -> v2``
    consuming ``vdiff`` of the block degree.

    Method (Molev--Sagan route): the verified top-block rule for
    ``prod_{j<=k}(x_j - t) G_{u1}`` gives a coefficient that is polynomial in the
    layer variable ``t``,

        ``C(t) = beta^(d - m) * prod_j f_j(t)``,   d = l(u2) - l(u1),

    with per-window-position factors ``f_j``: fixed value ``y`` ->
    ``(-y - t(1 + beta*y))/(1 + beta*y)``, left-moving persister -> ``1 + beta*t``,
    exit or right-mover -> ``1/(1 + beta*y)`` (``m`` = number of movers).  Since
    ``prod_{j<=k}(x_j - t) = sum_p e_p(x_1..x_k) (-t)^{k-p}``, extracting
    ``t``-coefficients of ``C`` gives the rule for each ``e_p``.  A v-path step of
    size ``vdiff`` is the z-side divided difference ``d_{vdiff} ... d_1``, which
    sends ``t^m`` to ``(-1)^vdiff h_{m - vdiff}`` in the ``vdiff + 1`` z-variables
    selected by ``call_zvars(v1, v2, k, i)`` (the same "screwed up" alphabet as
    the classical ``elem_sym_func``).  At ``vdiff = 0`` this is ``C(z_{v2(i)})``;
    at ``beta = 0`` it collapses to the classical ``elem_sym_func`` via
    ``E_{q,n}(y; z) = sum_j (-1)^j e_{q-j}(y) h_j(z)``.

    Returns a flat fraction ``(numer, {a: exp})`` with denominator
    ``prod_a (1 + beta*varl1[a])**exp``: the per-position factors are distributed
    so the numerator stays polynomial and the denominator is tracked as atom
    exponents, never as nested symbolic fractions.
    """
    from schubmult.symbolic.poly.schub_poly import call_zvars

    d = u2.inv - u1.inv
    window2 = [u2[j] for j in range(k)]
    coeffs = [S.One]
    denom = {}
    movers = 0
    for j in range(k):
        value = u1[j]
        if window2[j] == value:
            y = varl1[value]
            # (-y/(1+by) - t) = (-y - t(1+by)) / (1+by): polynomial numerator, atom bump.
            coeffs = _mul_linear(coeffs, -y, -(S.One + beta * y))
            denom[value] = denom.get(value, 0) + 1
        else:
            movers += 1
            if value in window2 and window2.index(value) < j:
                coeffs = _mul_linear(coeffs, S.One, beta)
            else:
                denom[value] = denom.get(value, 0) + 1
    if d < movers:
        return (S.Zero, {})
    zvars = [varl2[a] for a in call_zvars(v1, v2, k, i)][: vdiff + 1]
    total = S.Zero
    for m in range(vdiff, len(coeffs)):
        if coeffs[m] == S.Zero:
            continue
        total += coeffs[m] * _complete_homog(m - vdiff, zvars)
    sign = S.One if vdiff % 2 == 0 else -S.One
    return (sign * beta ** (d - movers) * total, denom)


def _groth_schub_vpath_mul(perm_dict, v, var2, var3, beta, as_frac=False):
    """``sum_u coeff_u G_u(x, var2) * S_v(x, var3)`` in the ``G`` basis.

    Mirrors ``schubmult_double``: expand ``S_v`` by ``compute_vpathdicts`` over the
    layers of ``theta(v^{-1})``, with ``elem_sym_perms`` -> marked-chain K-Pieri
    support and ``elem_sym_func`` -> ``_groth_elem_sym_frac``.

    All path sums are carried as flat fractions ``(numer, {atom: exp})`` and only
    reconstituted at the very end (or returned raw with ``as_frac=True``).
    """
    from schubmult.combinatorics.permutation import uncode
    from schubmult.utils.schub_lib import compute_vpathdicts

    v = Permutation(v)
    vn1 = ~v
    th = list(vn1.theta())
    while th and th[-1] == 0:
        th.pop()
    if not th:
        if as_frac:
            return {Permutation(w): (sympify(val), {}) for w, val in perm_dict.items()}
        return dict(perm_dict)
    mu = uncode(th)
    vmu = v * mu
    vpathdicts = compute_vpathdicts(tuple(th), vmu)
    ret_dict = {}
    for u, val in perm_dict.items():
        u = Permutation(u)
        vpathsums = {u: {Permutation([1, 2]): (sympify(val), {})}}
        for index in range(len(th)):
            k = th[index]
            newpathsums = {}
            for up, sums in vpathsums.items():
                for up2 in _top_block_support(up, k) | {up}:
                    for v_iter, steps in vpathdicts[index].items():
                        sumval = sums.get(v_iter)
                        if sumval is None or sumval[0] == S.Zero:
                            continue
                        for v2, vdiff, s in steps:
                            coeff = _groth_elem_sym_frac(k, index + 1, up, up2, v_iter, v2, vdiff, var2, var3, beta)
                            if coeff[0] == S.Zero:
                                continue
                            contrib = _frac_mul(sumval, coeff)
                            if s != 1:
                                contrib = (s * contrib[0], contrib[1])
                            bucket = newpathsums.setdefault(up2, {})
                            bucket[v2] = _frac_add(bucket.get(v2), contrib, var2, beta)
            vpathsums = newpathsums
        for ep, sums in vpathsums.items():
            pair = sums.get(vmu)
            if pair is not None and pair[0] != S.Zero:
                ret_dict[ep] = _frac_add(ret_dict.get(ep), pair, var2, beta)
    if as_frac:
        return {w: f for w, f in ret_dict.items() if f[0] != S.Zero}
    ret = {w: _frac_to_expr(f, var2, beta) for w, f in ret_dict.items()}
    return {w: coeff for w, coeff in ret.items() if coeff != S.Zero}


def grothmult_double(perm_dict, v, var2=None, var3=None, beta=None):
    r"""Multiply double Grothendieck polynomials, mirroring ``schubmult_double``.

    Computes the expansion of ``sum_u coeff_u G_u(x, var2) * G_v(x, var3)`` in
    the basis ``{G_w(x, var2)}`` and returns it as ``{w: coeff_w}``.

    ``v = s_k`` uses the verified chain formula of Corollary 8.2, and
    ``max_descent == 1`` folds that column by column.  General ``v`` goes through
    ``dgroth_to_dschub`` (exact, slow) and the conjectural vpath kernel
    ``_groth_schub_vpath_mul``, one run per double Schubert ``S_{v'}`` in the
    expansion of ``G_v``.

    The chain rank is inferred from the current permutation and selected positions.
    """
    if beta is None:
        beta = _default_beta
    var2 = _genset(var2)
    var3 = _genset(var3)

    v = Permutation(v)
    perm_dict = {Permutation(key): value for key, value in perm_dict.items()}
    if v.inv == 0:
        return perm_dict
    ret = {}
    for vprime, coeff in dgroth_to_dschub(v, var3, beta).items():
        for w, value in _groth_schub_vpath_mul(perm_dict, vprime, var2, var3, beta, as_frac=True).items():
            ret[w] = _frac_add(ret.get(w), (coeff * value[0], value[1]), var2, beta)
    out = {w: _frac_to_expr(f, var2, beta) for w, f in ret.items()}
    return {w: value for w, value in out.items() if value != S.Zero}
