"""Manifestly positive representations of double Schubert structure constants.

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
"""

from functools import cache

import psutil
from cachetools import cached
from cachetools.keys import hashkey

from schubmult.combinatorics.permutation import (
    Permutation,
    cycle,
    phi1,
    uncode,
)
from schubmult.symbolic import S, expand, prod, sympify, sympify_sympy, sympy_poly
from schubmult.symbolic.common_polys import _vars, efficient_subs, elem_sym_poly, schubpoly
from schubmult.utils.logging import get_logger
from schubmult.utils.schub_lib import (
    divdiffable,
    is_coeff_irreducible,
    is_split_two,
    pull_out_var,
    reduce_coeff,
    reduce_descents,
    try_reduce_u,
    try_reduce_v,
    will_formula_work,
)

from .double import schubmult_double, schubmult_double_pair, schubmult_double_pair_generic_alt

logger = get_logger(__name__)


def cbc_solver(msg=False):
    """PuLP ``COIN_CMD`` using the CBC binary bundled by ``pulp[cbc]``, falling back to ``cbc`` on PATH."""
    import pulp as pu

    try:
        from cbcbox import cbc_bin_path

        return pu.COIN_CMD(msg=msg, path=cbc_bin_path())
    except ImportError:
        return pu.COIN_CMD(msg=msg)


def compute_positive_rep(val, var2=None, var3=None, msg=False):
    """Express ``val`` as a nonnegative-integer combination of product-of-differences monomials.

    ``val`` must be a polynomial in ``var2``/``var3`` known (by positivity of
    double Schubert structure constants) to admit an expansion
    ``sum_b n_b * prod (var2_i - var3_j)`` with ``n_b >= 0`` integers. Builds a
    candidate spanning set of such product monomials from ``val``'s own
    monomials, then solves an integer program (via PuLP) for nonnegative
    integer coefficients ``n_b`` matching ``val`` exactly.

    Args:
        val: Symbolic polynomial expression in ``var2``/``var3``.
        var2: First secondary alphabet (``y``).
        var3: Second secondary alphabet (``z``).
        msg: Passed through to the LP solver as its ``msg`` (verbosity) option.

    Returns:
        A symbolic expression equal to ``val``, written as a sum of
        nonnegative-integer multiples of product-of-differences monomials.

    Raises:
        Exception: If the reconstructed expression does not equal ``val``
            (i.e. no valid nonnegative integer solution reproduces it exactly).
    """
    import pulp as pu

    try:
        return int(expand(val))
    except Exception:
        pass

    frees = val.free_symbols

    varsimp2 = [m for m in frees if var2.index(m) != -1]
    varsimp3 = [m for m in frees if var3.index(m) != -1]
    varsimp2.sort(key=lambda k: var2.index(k))
    varsimp3.sort(key=lambda k: var3.index(k))

    var22 = [sympify_sympy(v) for v in varsimp2]
    var33 = [sympify_sympy(v) for v in varsimp3]

    n1 = len(varsimp2)

    base_vectors = {}

    val_expr = expand(val)
    vec0 = {k: v for k, v in val_expr.subs({var3[1]: S.Zero}).as_coefficients_dict().items() if v != S.Zero}
    val_poly = sympy_poly(val_expr, *var22, *var33)

    mn = val_poly.monoms()
    L1 = tuple([0 for i in range(n1)])
    mn1L = []
    lookup = {}

    for mm0 in mn:
        key = mm0[n1:]
        if key not in lookup:
            lookup[key] = []
        mm0n1 = mm0[:n1]
        st = set(mm0n1)
        if len(st.intersection({0, 1})) == len(st) and 1 in st:
            lookup[key] += [mm0]
        if mm0n1 == L1:
            mn1L += [mm0]

    for mn1 in mn1L:
        comblistmn1 = [S.One]
        for i in range(n1, len(mn1)):
            if mn1[i] != 0:
                arr = [*comblistmn1]
                comblistmn12 = []
                mn1_2 = (*mn1[n1:i], 0, *mn1[i + 1 :])
                for mm0 in lookup[mn1_2]:
                    prd = sympify(
                        prod(
                            [varsimp2[k] - varsimp3[i - n1] for k in range(n1) if mm0[k] == 1],
                            start=S.One,
                        ),
                    )
                    comblistmn12 += [a * prd for a in arr]
                comblistmn1 = comblistmn12
        for i in range(len(comblistmn1)):
            b1 = comblistmn1[i]

            dct2 = {k: v for k, v in expand(b1).subs({var3[1]: S.Zero}).as_coefficients_dict().items() if v != S.Zero}
            bad = False
            for k in dct2:
                if abs(vec0.get(k, 0)) < abs(dct2[k]):
                    bad = True
                    break
            if not bad:
                base_vectors[b1] = dct2
    lp_prob = pu.LpProblem("Problem", pu.LpMinimize)
    vrs = {bv: lp_prob.add_variable(f"a{bv}", lowBound=0, cat="Integer") for bv in base_vectors}
    lp_prob += 0
    eqs = {}
    for bv, vec in base_vectors.items():
        for i in vec:
            bvi = int(vec[i])
            if bvi == 1:
                if i not in eqs:
                    eqs[i] = vrs[bv]
                else:
                    eqs[i] += vrs[bv]
            elif bvi != 0:
                if i not in eqs:
                    eqs[i] = bvi * vrs[bv]
                else:
                    eqs[i] += bvi * vrs[bv]
    for i in eqs:
        try:
            # PuLP >= 4 returns a bare False for `expr == <symengine Integer>`
            lp_prob += eqs[i] == int(vec0[i])
        except KeyError:
            raise

    try:
        solver = cbc_solver(msg)
        status = lp_prob.solve(solver)  # noqa: F841
    except KeyboardInterrupt:
        current_process = psutil.Process()
        children = current_process.children(recursive=True)
        for child in children:
            child_process = psutil.Process(child.pid)
            child_process.terminate()
            child_process.kill()
        raise KeyboardInterrupt()

    val2 = 0
    for k in base_vectors:
        x = vrs[k].value()
        # round, don't truncate: solvers return near-integers like 0.9999999999996
        if x is not None and round(x) != 0:
            val2 += round(x) * k
    if expand(val - val2, func=True) != 0:
        raise ValueError("Failed to find a positive representation.")
    return val2


@cached(
    cache={},
    key=lambda val, u2, v2, w2, var2=None, var3=None, msg=False, sign_only=False, optimize=True: hashkey(val, u2, v2, w2, var2, var3, msg, sign_only, optimize),
)
def posify(
    val,
    u2,
    v2,
    w2,
    var2=None,
    var3=None,
    msg=False,
    sign_only=False,
    optimize=True,
    n=_vars.n,
):
    """Manifestly positive representation of the structure constant ``c^{w2}_{u2,v2}(var2, var3)``.

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

    Args:
        val: The structure constant to re-express positively.
        u2: First factor's permutation.
        v2: Second factor's permutation.
        w2: Target permutation (coefficient of ``S_{w2}``).
        var2: First secondary alphabet.
        var3: Second secondary alphabet.
        msg: Verbosity flag passed down to ``compute_positive_rep``'s LP solver.
        sign_only: If ``True``, only determine and return the sign of ``val``
            (``-1``, ``0``, or ``1``) rather than a full positive expression.
        optimize: If ``False``, return ``val`` unchanged when no closed-form
            reduction applies (skip the LP fallback); if ``None`` and that
            case is reached, raise.
        n: Size of the ambient alphabet used when no other bound is available.

    Returns:
        A manifestly positive expression equal to ``val`` (or, if
        ``sign_only``, one of ``-1``, ``0``, ``1``).
    """
    if not v2.has_pattern([1, 4, 2, 3]) and not v2.has_pattern([4, 1, 3, 2]) and not v2.has_pattern([3, 1, 4, 2]) and not v2.has_pattern([1, 4, 3, 2]):
        logger.debug("Recording new characterization was used")
        return schubmult_double({u2: 1}, v2, var2, var3).get(w2, 0)
    oldval = val
    if u2.inv + v2.inv - w2.inv == 0:
        return val

    if not sign_only and expand(val) == 0:
        return 0

    u, v, w = u2, v2, w2
    if is_coeff_irreducible(u2, v2, w2):
        u, v, w = try_reduce_u(u2, v2, w2)
        if is_coeff_irreducible(u, v, w):
            u, v, w = u2, v2, w2
            if is_coeff_irreducible(u, v, w):
                w0 = w
                u, v, w = reduce_descents(u, v, w)
                if is_coeff_irreducible(u, v, w):
                    u, v, w = reduce_coeff(u, v, w)
                    if is_coeff_irreducible(u, v, w):
                        while is_coeff_irreducible(u, v, w) and w0 != w:
                            w0 = w
                            u, v, w = reduce_descents(u, v, w)
                            if is_coeff_irreducible(u, v, w):
                                u, v, w = reduce_coeff(u, v, w)

    if w != w2 and sign_only:
        return 0

    if is_coeff_irreducible(u, v, w):
        u3, v3, w3 = try_reduce_v(u, v, w)
        if not is_coeff_irreducible(u3, v3, w3):
            u, v, w = u3, v3, w3
        else:
            u3, v3, w3 = try_reduce_u(u, v, w)
            if not is_coeff_irreducible(u3, v3, w3):
                u, v, w = u3, v3, w3
    _split_two_b, _split_two = is_split_two(u, v, w)

    if len([i for i in v.code if i != 0]) == 1:
        if sign_only:
            return 0
        cv = v.code
        for i in range(len(cv)):
            if cv[i] != 0:
                k = i + 1
                p = cv[i]
                break
        inv_u = u.inv
        r = w.inv - inv_u
        val = 0
        w2 = w
        hvarset = [w2[i] for i in range(min(len(w2), k))] + [i + 1 for i in range(len(w2), k)] + [w2[b] for b in range(k, len(u)) if u[b] != w2[b]] + [w2[b] for b in range(len(u), len(w2))]

        return elem_sym_poly(
            p - r,
            k + p - 1,
            [-var3[i] for i in range(1, n)],
            [-var2[i] for i in hvarset],
        )

    if will_formula_work(v, u) or u.dominates(w):
        if sign_only:
            return 0
        return dualcoeff(u, v, w, var2, var3)

    if not v.has_pattern([1, 4, 2, 3]) and not v.has_pattern([4, 1, 3, 2]) and not v.has_pattern([3, 1, 4, 2]) and not v.has_pattern([1, 4, 3, 2]):
        logger.debug("Recording new characterization was used")
        return schubmult_double({u: 1}, v, var2, var3).get(w, 0)

    if w.inv - u.inv == 1:
        if sign_only:
            return 0
        a, b = -1, -1
        for i in range(len(w)):
            if a == -1 and u[i] != w[i]:
                a = i
            elif (i >= len(u) and w[i] != i + 1) or (b == -1 and u[i] != w[i]):
                b = i
        arr = [[[], v]]
        d = -1
        for i in range(len(v) - 1):
            if v[i] > v[i + 1]:
                d = i + 1
        for i in range(d):
            arr2 = []
            if i in [a, b]:
                continue
            i2 = 1
            if i > b:
                i2 += 2
            elif i > a:
                i2 += 1
            for vr, v2 in arr:
                dpret = pull_out_var(i2, v2)
                for v3r, v3 in dpret:
                    arr2 += [[[*vr, v3r], v3]]
            arr = arr2
        val = 0
        for L in arr:
            v3 = L[-1]
            if v3[0] < v3[1]:
                continue
            v3 = v3.swap(0, 1)
            toadd = 1
            for i in range(d):
                if i in [a, b]:
                    continue
                i2 = i
                if i > b:
                    i2 = i - 2
                elif i > a:
                    i2 = i - 1
                oaf = L[0][i2]
                if i >= len(w):
                    yv = i + 1
                else:
                    yv = w[i]
                for j in range(len(oaf)):
                    toadd *= var2[yv] - var3[oaf[j]]
            toadd *= schubpoly(v3, [0, var2[w[a]], var2[w[b]]], var3)
            val += toadd
        return val

    if will_formula_work(u, v):
        if sign_only:
            return 0
        return forwardcoeff(u, v, w, var2, var3)

    c1 = (~u).code
    c2 = (~w).code

    if u.one_dominates(w):
        if sign_only:
            return 0
        while c1[0] != c2[0]:
            w = w.swap(c2[0] - 1, c2[0])
            v = v.swap(c2[0] - 1, c2[0])

            c2 = (~w).code

        if c1[0] == c2[0]:
            if sign_only:
                return 0
            vp = pull_out_var(c1[0] + 1, v)
            u3 = phi1(u)
            w3 = phi1(w)
            val = 0
            for arr, v3 in vp:
                tomul = prod([var2[1] - var3[arr[i]] for i in range(len(arr))])

                val2 = schubmult_double_pair(u3, v3, var2, var3).get(
                    w3,
                    0,
                )
                val2 = posify(val2, u3, v3, w3, var2[1:], var3, msg, optimize=optimize)
                val += tomul * val2

            return val

    if not sign_only:
        if optimize:
            if u.inv + v.inv - w.inv == 1:
                val2 = compute_positive_rep(val, var2, var3, msg)
            else:
                val2 = compute_positive_rep(val, var2, var3, msg)
            if val2 is not None:
                val = val2
            return val
        if optimize is None:
            raise ValueError("Optimize is None but no positive formula applied")
        return oldval
    d = expand(val).as_coefficients_dict()
    for v in d.values():
        if v < 0:
            return -1
    return 1


def shiftsub(pol, var2=None):
    """Shift every ``var2[i]`` in ``pol`` up to ``var2[i + 1]`` (for ``i`` in ``0..98``)."""
    subs_dict = {var2[i]: var2[i + 1] for i in range(99)}
    return efficient_subs(sympify(pol), subs_dict)


def posify_generic_partial(val, u2, v2, w2):
    """``posify`` specialized to the fixed generic alphabets ``_vars.var_g1``/``_vars.var_g2``.

    Asserts (raises on mismatch) that the recomputed positive expression equals
    the input ``val``, as a consistency check.
    """
    val2 = val
    val = posify(val, u2, v2, w2, var2=_vars.var_g1, var3=_vars.var_g2, msg=True, sign_only=False, optimize=False)
    if expand(val - val2) != 0:
        raise Exception(f"{val=} {val2=} {u2=} {v2=} {w2=}")

    return val


@cache
def schubmult_generic_partial_posify(u2, v2):
    """Manifestly positive expansion of ``S_{u2}(x, var_g1) * S_{v2}(x, var_g2)``.

    Returns ``{w2: coeff}`` where each ``coeff`` is the positive representation
    (via ``posify_generic_partial``) of the corresponding
    ``schubmult_double_pair_generic_alt`` coefficient.
    """
    return {w2: posify_generic_partial(val, u2, v2, w2) for w2, val in schubmult_double_pair_generic_alt(u2, v2).items()}


def forwardcoeff(u, v, perm, var2=None, var3=None):
    """Closed-form structure constant ``c^{perm}_{u,v}(var2, var3)`` for the "forward" case
    (used when ``will_formula_work(u, v)`` holds in ``posify``).

    Writes ``muv = uncode(v.theta())`` and reduces to a lookup in
    ``schubmult_double_pair(u, muv, var2, var3)`` when the length condition
    ``(perm * (~v * muv)).inv == (~v * muv).inv + perm.inv`` holds; returns 0 otherwise.
    """
    th = v.theta()
    muv = uncode(th)
    vmun1 = (~v) * muv

    w = perm * vmun1
    if w.inv == vmun1.inv + perm.inv:
        coeff_dict = schubmult_double_pair(u, muv, var2, var3)
        return coeff_dict.get(w, 0)
    return 0


def dualcoeff(u, v, perm, var2=None, var3=None):
    """Closed-form structure constant ``c^{perm}_{u,v}(var2, var3)`` for the "dual" case
    (used in ``posify`` when ``will_formula_work(v, u)`` holds or ``u`` dominates ``perm``).

    When ``u`` is the identity, reduces directly to a single Schubert
    polynomial ``schubpoly(v * (~perm), var2, var3)``. Otherwise expands via
    ``dualpieri`` (directly if ``u`` dominates ``perm``, or after rewriting to
    ``u``'s dominant permutation ``uncode(u.theta())`` otherwise), summing
    products of ``(var2_{i+1} - var3_j)`` factors against a final
    ``schubpoly`` term.
    """
    if u.inv == 0:
        vp = v * (~perm)
        if vp.inv == v.inv - perm.inv:
            return schubpoly(vp, var2, var3)
    dpret = []
    ret = 0
    if u.dominates(perm):
        dpret = dualpieri(u, v, perm)
    else:
        dpret = []
        th = u.theta()
        muu = uncode(th)
        umun1 = (~u) * muu
        w = perm * umun1

        if w.inv == umun1.inv + perm.inv:
            dpret = dualpieri(muu, v, w)

    for vlist, vp in dpret:
        toadd = 1
        for i in range(len(vlist)):
            for j in range(len(vlist[i])):
                toadd *= var2[i + 1] - var3[vlist[i][j]]
        toadd *= schubpoly(vp, var2, var3, len(vlist) + 1)
        ret += toadd
    return ret


def dualpieri(mu, v, w):
    """Dual Pieri expansion used by ``dualcoeff``: enumerate the data witnessing
    ``S_mu * S_v -> S_w`` when ``mu`` is dominant.

    Compares ``mu``'s inverse code against ``w``'s inverse code layer by layer,
    peeling one "cycle" of variables per layer via ``divdiffable``/``pull_out_var``,
    and returns the list of ``[vlist, vp]`` pairs consumed by ``dualcoeff`` to
    build the final positive expression (empty list if ``w`` is not reachable
    from ``mu``, ``v`` this way).
    """
    lm = (~mu).code
    cn1w = (~w).code
    while len(lm) > 0 and lm[-1] == 0:
        lm.pop()
    while len(cn1w) > 0 and cn1w[-1] == 0:
        cn1w.pop()
    if len(cn1w) < len(lm):
        return []
    for i in range(len(lm)):
        if lm[i] > cn1w[i]:
            return []
    c = Permutation([1, 2])
    for i in range(len(lm), len(cn1w)):
        c = cycle(i - len(lm) + 1, cn1w[i]) * c

    res = [[[], v]]

    for i in range(len(lm)):
        res2 = []
        for vlist, vplist in res:
            vp = vplist
            vpl = divdiffable(vp, cycle(lm[i] + 1, cn1w[i] - lm[i]))
            if len(vpl) == 0:
                continue
            vl = pull_out_var(lm[i] + 1, vpl)
            for pw, vpl2 in vl:
                res2 += [[[*vlist, pw], vpl2]]
        res = res2
    if len(lm) == len(cn1w):
        return res
    res2 = []
    for vlist, vplist in res:
        vp = vplist
        vpl = divdiffable(vp, c)
        if len(vpl) == 0:
            continue
        res2 += [[vlist, vpl]]
    return res2
