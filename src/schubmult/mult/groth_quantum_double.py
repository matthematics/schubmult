r"""Quantum double Grothendieck multiplication: the Molev--Sagan machinery for ``QK_T(Fl_n)``.

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
"""

from functools import cache
from itertools import combinations, product
from math import comb

from schubmult.abc import beta as _default_beta
from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.mult.groth_double import (
    _frac_add,
    _frac_mul,
    _frac_to_expr,
    _genset,
    _groth_elem_sym_frac,
    _rank,
    _top_block_coeff,
    dgroth_to_dschub,
    groth_elem_sym_poly,
)
from schubmult.symbolic import S, prod, sympify, sympify_sympy
from schubmult.symbolic.poly.schub_poly import _vars
from schubmult.utils.perm_utils import add_perm_dict
from schubmult.utils.schub_lib import compute_vpathdicts

__all__ = [
    "groth_elem_sym_poly_q",
    "grothmult_q_double",
    "grothmult_q_double_dict",
    "grothmult_q_double_pieri",
    "grothmult_q_double_top",
    "lm_quantize",
    "qgroth_poly",
    "quantum_elem_sym",
    "quantum_pieri_chains",
]


def _ns_prec(label1, label2):
    """Naito--Sagaki order on labels: ``(a, b) < (c, d)`` iff ``b > d``, or ``b == d`` and ``a < c``."""
    (a, b), (c, d) = label1, label2
    return b > d or (b == d and a < c)


@cache
def quantum_pieri_chains(u, k):
    r"""Endpoints of the Naito--Sagaki ``k``-Pieri chains from ``u`` in the quantum Bruhat graph.

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
    """
    u = Permutation(u)
    rank = _rank(u, tuple(range(1, k + 1)))
    labels = [(a, b) for b in range(rank, k, -1) for a in range(1, k + 1)]
    found = {}

    def walk(w, path, dvec):
        found.setdefault(w, set()).add((len(path), dvec))
        repeated_last = bool(path) and any(t[0] == path[-1][0] for t in path[:-1])
        for label in labels:
            a, b = label
            if label in path or (path and b > path[-1][1]):
                continue
            if repeated_last and not _ns_prec(path[-1], label):
                continue
            stepped = w.swap(a - 1, b - 1)
            if stepped.inv == w.inv + 1:
                new_d = dvec
            elif stepped.inv == w.inv - (2 * (b - a) - 1):
                new_d = tuple(dvec[j] + (1 if a <= j + 1 < b else 0) for j in range(rank - 1))
            else:
                continue
            walk(stepped, [*path, label], new_d)

    walk(u, [], (0,) * (rank - 1))
    out = {}
    for w, data in found.items():
        if len(data) != 1:
            raise ValueError(f"ambiguous quantum chain data for u={list(u)}, k={k}, w={list(w)}: {sorted(data)}")
        out[w] = next(iter(data))
    return out


@cache
def _qmon(dvec, q_var):
    """``prod q_j ** D[j - 1]`` for an exponent tuple ``D`` (memoized)."""
    return prod([q_var[j + 1] ** e for j, e in enumerate(dvec) if e])


def _fate(u, w, k):
    """Sort the window values ``u(i)``, ``i <= k``, by their fate in ``w``: ``(fixed, left, out)`` value lists."""
    window = [w[i] for i in range(k)]
    fixed, left, out = [], [], []
    for i in range(k):
        value = u[i]
        if window[i] == value:
            fixed.append(value)
        elif value in window and window.index(value) < i:
            left.append(value)
        else:
            out.append(value)
    return fixed, left, out


def grothmult_q_double_top(coeff_dict, k, zvar=None, var2=None, beta=None, q_var=None):
    r"""Multiply ``sum_u coeff_u G^q_u(x, var2)`` by the quantized top block ``Q(prod_{i<=k}(x_i + zvar))``.

    The multiplier is ``groth_elem_sym_poly_q(k, k, zvar, x, beta, q_var, fgl=False)``, i.e.
    ``beta^{-k} sum_l (-(1 - beta z))^{k-l} F^k_l(X)``.  The coefficient of ``G^q_w`` is
    ``_top_block_coeff`` evaluated with the quantum chain length, times the quantum weight
    ``q^D``; see the module docstring.
    """
    if beta is None:
        beta = _default_beta
    if zvar is None:
        zvar = S.Zero
    if q_var is None:
        q_var = _vars.q_var
    var2 = _genset(var2)

    ret = {}
    for u, val in coeff_dict.items():
        u = Permutation(u)
        if k == 0:
            ret[u] = ret.get(u, S.Zero) + val
            continue
        for w, (length, dvec) in quantum_pieri_chains(u, k).items():
            ret[w] = ret.get(w, S.Zero) + val * _qmon(dvec, q_var) * _top_block_coeff(u, w, k, zvar, var2, beta, length=length)
    return ret


def groth_elem_sym_poly_q(p, k, zvar, var_x, beta, q_var=None, fgl=True):
    r"""Quantization of ``groth_elem_sym_poly``: ``Q(e_p(x_1 (+) z, ..., x_k (+) z))``.

    ``fgl=False`` quantizes the plain ``e_p(x_1 + z, ..., x_k + z)`` instead, whose ``p = k``
    case is the top block of ``grothmult_q_double_top``.  Both are symmetric in ``x_1..x_k`` of
    degree ``<= 1`` in each variable, so ``lm_quantize`` with ``k + 1`` slots applies.
    """
    if q_var is None:
        q_var = _vars.q_var
    if fgl:
        poly = groth_elem_sym_poly(p, k, zvar, var_x, beta)
    else:
        acc = [S.One] + [S.Zero] * p
        for i in range(1, k + 1):
            for r in range(p, 0, -1):
                acc[r] = acc[r] + acc[r - 1] * (var_x[i] + zvar)
        poly = acc[p]
    return lm_quantize(poly, k + 1, var_x, beta, q_var)


def _q_pieri_coeff(u, w, k, p, zvar, var2, beta, length, fgl):
    r"""Coefficient of ``G^q_w`` in ``Q(e_p(x (+) z)) G^q_u`` (without the ``q^D`` factor).

    With ``T(t) = prod_{Fix}(t + (-)y_a) (1 - beta t)^{|Left|} = sum_r T_r t^r`` the top block
    gives ``e_j(x) G^q_u -> beta^{len - m} prod_{Out}(1 + beta y)^{-1} T_{k-j} G^q_w``, and
    ``e_p(x (+) z) = sum_j binom(k-j, p-j) z^{p-j} (1 + beta z)^j e_j(x)``.
    """
    fixed, left, out = _fate(u, w, k)
    m = len(left) + len(out)
    if length < m:
        return S.Zero
    poly = [S.One]
    for a in fixed:
        root = -var2[a] / (S.One + beta * var2[a])
        poly = [(poly[r - 1] if r else S.Zero) + root * (poly[r] if r < len(poly) else S.Zero) for r in range(len(poly) + 1)]
    for _ in left:
        poly = [(poly[r] if r < len(poly) else S.Zero) - beta * (poly[r - 1] if r else S.Zero) for r in range(len(poly) + 1)]
    scale = S.One + beta * zvar if fgl else S.One
    total = S.Zero
    for j in range(len(out), p + 1):
        total += comb(k - j, p - j) * zvar ** (p - j) * scale**j * poly[k - j]
    if total == S.Zero:
        return S.Zero
    return beta ** (length - m) * prod([S.One / (S.One + beta * var2[a]) for a in out]) * total


def grothmult_q_double_pieri(coeff_dict, p, k, zvar=None, var2=None, beta=None, q_var=None, fgl=True):
    r"""Multiply ``sum_u coeff_u G^q_u(x, var2)`` by ``groth_elem_sym_poly_q(p, k, zvar, x, beta, q_var, fgl)``.

    Closed form from the top block (see ``_q_pieri_coeff``); ``fgl=False`` multiplies by the
    quantization of ``e_p(x_1 + z, ..., x_k + z)``.  At ``q = 0`` this agrees with the exact
    fold ``grothmult_double_pieri``.
    """
    if beta is None:
        beta = _default_beta
    if zvar is None:
        zvar = S.Zero
    if q_var is None:
        q_var = _vars.q_var
    var2 = _genset(var2)

    ret = {}
    for u, val in coeff_dict.items():
        u = Permutation(u)
        if p == 0:
            ret[u] = ret.get(u, S.Zero) + val
            continue
        for w, (length, dvec) in quantum_pieri_chains(u, k).items():
            coeff = _q_pieri_coeff(u, w, k, p, zvar, var2, beta, length, fgl)
            if coeff != S.Zero:
                ret[w] = ret.get(w, S.Zero) + val * _qmon(dvec, q_var) * coeff
    return ret


def _qgroth_schub_vpath_mul(perm_dict, v, var2, var3, beta, q_var, as_frac=False):
    """``sum_u coeff_u G^q_u(x, var2) * S^q_v(x, var3)`` in the ``G^q`` basis.

    Mirrors ``schubmult_q_double`` (strict theta, one layer per level) with the quantum Pieri
    support ``quantum_pieri_chains`` and the coefficient ``_groth_elem_sym_frac`` at the quantum
    chain length, times ``q^D``.  Path sums are flat fractions ``(numer, {atom: exp})`` as in
    ``_groth_schub_vpath_mul``.
    """
    v = Permutation(v)
    th = list((~v).strict_theta())
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
        for index, k in enumerate(th):
            layer = vpathdicts[index]
            i = index + 1
            newpathsums = {}
            for up, sums in vpathsums.items():
                live = [(v_iter, sumval, layer[v_iter]) for v_iter, sumval in sums.items() if sumval[0] != S.Zero and v_iter in layer]
                if not live:
                    continue
                for up2, (length, dvec) in quantum_pieri_chains(up, k).items():
                    qmon = _qmon(dvec, q_var)
                    bucket = None
                    for v_iter, sumval, steps in live:
                        for v2, vdiff, s in steps:
                            coeff = _groth_elem_sym_frac(k, i, up, up2, v_iter, v2, vdiff, var2, var3, beta, length=length)
                            if coeff[0] == S.Zero:
                                continue
                            contrib = _frac_mul(sumval, (s * qmon * coeff[0], coeff[1]))
                            if bucket is None:
                                bucket = newpathsums.setdefault(up2, {})
                            bucket[v2] = _frac_add(bucket.get(v2), contrib, var2, beta)
            vpathsums = newpathsums
        for ep, sums in vpathsums.items():
            pair = sums.get(vmu)
            if pair is not None and pair[0] != S.Zero:
                ret_dict[ep] = _frac_add(ret_dict.get(ep), pair, var2, beta)
    if as_frac:
        return {w: f for w, f in ret_dict.items() if f[0] != S.Zero}
    out = {w: _frac_to_expr(f, var2, beta) for w, f in ret_dict.items()}
    return {w: coeff for w, coeff in out.items() if coeff != S.Zero}


def grothmult_q_double(perm_dict, v, var2=None, var3=None, beta=None, q_var=None):
    r"""Multiply quantum double Grothendieck polynomials, mirroring ``schubmult_q_double``.

    Returns the expansion of ``sum_u coeff_u G^q_u(x, var2) * G^q_v(x, var3)`` in the basis
    ``{G^q_w(x, var2)}`` as ``{w: coeff_w}``, coefficients rational in ``var2`` (denominators
    are products of ``1 + beta*var2[a]``) and polynomial in ``var3``, ``beta``, ``q``.

    ``G^q_v(x, var3)`` is expanded through ``dgroth_to_dschub`` and one run of the quantum
    v-path kernel ``_qgroth_schub_vpath_mul`` per double Schubert term.  At ``q = 0`` this is
    ``grothmult_double``; at ``beta = 0`` it is ``schubmult_q_double``.
    """
    if beta is None:
        beta = _default_beta
    if q_var is None:
        q_var = _vars.q_var
    var2 = _genset(var2)
    var3 = _genset(var3)

    v = Permutation(v)
    perm_dict = {Permutation(key): value for key, value in perm_dict.items()}
    if v.inv == 0:
        return perm_dict
    ret = {}
    for vprime, coeff in dgroth_to_dschub(v, var3, beta).items():
        for w, value in _qgroth_schub_vpath_mul(perm_dict, vprime, var2, var3, beta, q_var, as_frac=True).items():
            ret[w] = _frac_add(ret.get(w), (coeff * value[0], value[1]), var2, beta)
    out = {w: _frac_to_expr(f, var2, beta) for w, f in ret.items()}
    return {w: coeff for w, coeff in out.items() if coeff != S.Zero}


def grothmult_q_double_dict(perm_dict1, perm_dict2, var2=None, var3=None, beta=None, q_var=None):
    """Product of two coefficient dicts: ``sum_v coeff2_v grothmult_q_double(perm_dict1, v, ...)``."""
    ret = {}
    for v, coeff in perm_dict2.items():
        ret = add_perm_dict(ret, {w: coeff * value for w, value in grothmult_q_double(perm_dict1, v, var2, var3, beta, q_var).items()})
    return ret


def quantum_elem_sym(l, k, var_x, beta, q_var=None):
    r"""``F^k_l(X) = sum_{J in [k], |J| = l} prod_{j in J, j+1 not in J} (1 - beta^2 q_j) prod_{j in J} X_j``, ``X_j = 1 + beta*x_j``.

    The Lenart--Maeno quantization of ``e_l(X_1..X_k)``; ``beta = -1`` gives the ``F^k_l`` of
    Maeno--Naito--Sagaki with ``Q_j = q_j``.  No ``1 - Q_N := 1`` convention is applied: that
    belongs to the defining ideal of ``QK_T(Fl_N)``, not to the polynomials.
    """
    if q_var is None:
        q_var = _vars.q_var
    total = S.Zero
    for J in combinations(range(1, k + 1), l):
        Js = set(J)
        term = S.One
        for j in J:
            term *= S.One + beta * var_x[j]
            if j + 1 not in Js:
                term *= S.One - beta**2 * q_var[j]
        total += term
    return total


@cache
def _sem_basis(N):
    """Standard elementary monomial exponent tuples ``(i_1..i_{N-1})``, ``i_j <= j``, and the
    inverse (as a rational ``sympy.Matrix``) of the matrix expressing them in the monomials
    ``X^a``, ``a_i <= N - i`` (cached).
    """
    import sympy
    from sympy.polys.matrices import DomainMatrix

    X = sympy.symbols(f"X1:{N}")
    monos = sorted(product(*[range(N - i + 1) for i in range(1, N)]))
    mono_idx = {m: r for r, m in enumerate(monos)}
    sems = sorted(product(*[range(j + 1) for j in range(1, N)]))
    rows = [[sympy.QQ.zero] * len(sems) for _ in monos]
    for c, sem in enumerate(sems):
        expr = sympy.Integer(1)
        for j, i_j in enumerate(sem, start=1):
            expr *= sum(sympy.prod(J) for J in combinations(X[:j], i_j))
        for m, coeff in sympy.Poly(sympy.expand(expr), *X).terms():
            rows[mono_idx[m]][c] = sympy.QQ.convert(coeff)
    inv = DomainMatrix(rows, (len(monos), len(sems)), sympy.QQ).inv().to_Matrix().tolist()
    return X, monos, sems, inv


def lm_quantize(poly, N, var_x, beta, q_var=None):
    r"""Lenart--Maeno quantization of a polynomial in ``x_1..x_{N-1}`` of degree ``<= N - i`` in ``x_i``.

    Expands ``poly`` in the standard elementary monomials ``prod_j e_{i_j}(X_1..X_j)`` of
    ``X_i = 1 + beta*x_i`` (a basis of that span; ``_sem_basis``) and replaces each
    ``e_{i_j}(X_1..X_j)`` by ``quantum_elem_sym(i_j, j)``.  Linear over everything but ``x``,
    stable in ``N`` (``poly`` may use fewer variables), and ``Q_j = beta^2 q_j`` so the result
    is polynomial in ``beta``.  Raises if ``poly`` is not in the span.
    """
    import sympy

    if q_var is None:
        q_var = _vars.q_var
    X, monos, sems, inv = _sem_basis(N)
    bs = sympify_sympy(beta)
    sub = {sympify_sympy(var_x[i]): (X[i - 1] - 1) / bs for i in range(1, N)}
    f = sympy.expand(sympify_sympy(poly).xreplace(sub))
    P = sympy.Poly(f, *X)
    coeffs = [sympy.Integer(0)] * len(monos)
    mono_idx = {m: r for r, m in enumerate(monos)}
    for m, c in P.terms():
        if m not in mono_idx:
            raise ValueError(f"monomial X^{m} outside the standard elementary monomial span for N={N}")
        coeffs[mono_idx[m]] = c
    total = S.Zero
    for r, sem in enumerate(sems):
        c = sum((inv[r][s] * coeffs[s] for s in range(len(monos)) if coeffs[s] != 0), sympy.Integer(0))
        if c == 0:
            continue
        term = sympify(sympy.cancel(c))
        for j, i_j in enumerate(sem, start=1):
            if i_j:
                term *= quantum_elem_sym(i_j, j, var_x, beta, q_var)
        total += term
    return sympify(sympy.expand(sympy.cancel(sympify_sympy(total))))


def qgroth_poly(v, var_x=None, var_y=None, beta=None, q_var=None):
    r"""The quantum double Grothendieck polynomial ``G^q_v(var_x; var_y)`` as an explicit expression.

    ``lm_quantize`` applied to ``grothendieck_poly(v)`` with ``N = len(v)`` slots.  At
    ``beta = -1`` this is the Maeno--Naito--Sagaki ``G^Q_v(x, y)`` with ``Q_j = q_j``.  Slow
    (symbolic); intended for verification.
    """
    from schubmult.symbolic.poly.schub_poly import grothendieck_poly
    from schubmult.symbolic.poly.variables import GeneratingSet

    if beta is None:
        beta = _default_beta
    if q_var is None:
        q_var = _vars.q_var
    var_x = GeneratingSet("x") if var_x is None else _genset(var_x)
    var_y = GeneratingSet("y") if var_y is None else _genset(var_y)
    v = Permutation(v)
    return lm_quantize(grothendieck_poly(v, var_x, var_y, beta), max(len(v), 2), var_x, beta, q_var)
