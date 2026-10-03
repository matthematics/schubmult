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
from schubmult.symbolic import S, expand, prod, sympify
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


def _monomial_dict(expr, variables):
    """``{exponent tuple: int}`` of the expanded polynomial ``expr`` in ``variables`` (SymEngine, no SymPy)."""
    pos = {v: i for i, v in enumerate(variables)}
    out = {}
    for mono, c in expand(expr).as_coefficients_dict().items():
        e = [0] * len(variables)
        if mono != 1:
            for a in mono.args if mono.is_Mul else (mono,):
                if a.is_Pow:
                    base, k = a.args
                    e[pos[base]] += int(k)
                else:
                    e[pos[a]] += 1
        out[tuple(e)] = out.get(tuple(e), 0) + int(c)
    return {m: c for m, c in out.items() if c}


def _vec_mul(a, b):
    """Product of two sparse polynomials in exponent-tuple form."""
    out = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = tuple(x + y for x, y in zip(m1, m2))
            out[m] = out.get(m, 0) + c1 * c2
    return {m: c for m, c in out.items() if c}


def _difference_vec(ks, j, n1, n3):
    """``prod_{k in ks} (y_k + z_j)`` in exponent-tuple form (``y`` positions ``0..n1-1``, ``z_j`` at ``n1 + j``).

    All coefficients are positive: this is the product of differences ``prod (y_k - z_j)`` with the
    sign ``(-1)^(z-degree)`` of each monomial stripped, the convention used throughout the solver.
    """
    vec = {(0,) * (n1 + n3): 1}
    for k in ks:
        new = {}
        for m, c in vec.items():
            m1 = list(m)
            m1[k] += 1
            new[tuple(m1)] = new.get(tuple(m1), 0) + c
            m2 = list(m)
            m2[n1 + j] += 1
            new[tuple(m2)] = new.get(tuple(m2), 0) + c
        vec = new
    return vec


def _fits(vec, target):
    return all(target.get(m, 0) >= c for m, c in vec.items())


def _enumerate_candidates(target, n1, n3):
    """Products of differences that can occur in a positive decomposition of ``target`` (homogeneous,
    sign-stripped, nonnegative), as ``(pairs, vec)``: the factors ``(k, j)`` and the positive vector.

    A candidate is a bipartite multigraph between ``y`` and ``z`` whose pure-``z`` monomial ``z^beta``
    is a monomial of ``target``; its ``z_j`` group is a subset ``S`` of the ``y``'s with ``|S| = beta_j``
    and ``y^S z^(beta - beta_j e_j)`` a monomial of ``target``. Groups are assigned one ``z`` at a time,
    and a partial product is discarded as soon as one of its monomials (times the pure-``z`` part of the
    groups still to come, a monomial every completion has with at least that coefficient) exceeds the
    target: no complete candidate below it can occur.
    """
    by_z = {}
    for m in target:
        by_z.setdefault(m[n1:], []).append(m[:n1])
    zero = (0,) * (n1 + n3)
    cands = []
    for beta in by_z:
        if zero[:n1] not in by_z[beta]:  # z^beta itself must be a monomial of the target
            continue
        groups = [j for j in range(n3) if beta[j]]
        options = []
        for j in groups:
            base = tuple(0 if i == j else beta[i] for i in range(n3))
            options.append([tuple(k for k in range(n1) if y[k]) for y in by_z.get(base, []) if set(y) <= {0, 1} and sum(y) == beta[j]])
        if not all(options):
            continue
        suffix = []  # pure-z exponent tuple of the groups after each position
        for gi in range(len(groups)):
            s = [0] * (n1 + n3)
            for j in groups[gi + 1 :]:
                s[n1 + j] = beta[j]
            suffix.append(tuple(s))

        def rec(gi, pairs, vec):
            if gi == len(groups):
                cands.append((pairs, vec))
                return
            j, suf = groups[gi], suffix[gi]
            for ks in options[gi]:
                new = _vec_mul(vec, _difference_vec(ks, j, n1, n3))
                if all(target.get(tuple(a + b for a, b in zip(m, suf)), 0) >= c for m, c in new.items()):
                    rec(gi + 1, pairs + tuple((k, j) for k in ks), new)

        rec(0, (), {zero: 1})
    return cands


def _peel(vecs, resid, alive):
    """Subtract the forced candidates in place: while some monomial of ``resid`` is covered by exactly
    one candidate (among ``alive``) that still fits, that candidate's multiplicity is determined.
    Returns ``(forced {index: n}, surviving indices)``; ``None`` if a monomial of the residual is
    covered by no candidate (no decomposition over ``vecs``).
    """
    forced = {}
    while True:
        alive = [i for i in alive if _fits(vecs[i], resid)]
        cover = {}
        for i in alive:
            for m in vecs[i]:
                cover.setdefault(m, []).append(i)
        progress = False
        for m, c in resid.items():
            if c == 0:
                continue
            idx = cover.get(m)
            if not idx:
                return None
            if len(idx) == 1:
                i = idx[0]
                vec = vecs[i]
                if c % vec[m]:
                    return None
                n = c // vec[m]
                if not all(resid.get(mm, 0) >= n * cc for mm, cc in vec.items()):
                    return None
                for mm, cc in vec.items():
                    resid[mm] -= n * cc
                forced[i] = forced.get(i, 0) + n
                progress = True
                break
        if not progress:
            for m in [m for m, c in resid.items() if c == 0]:
                del resid[m]
            return forced, alive


def _solve_positive_system(vecs, target, msg):
    """Nonnegative integers ``x`` with ``sum_i x_i vecs[i] == target`` (all entries nonnegative), or ``None``.

    First the exact search (:func:`_cover_search`) within a node budget. If that runs out, rounds of:
    solve the LP relaxation and fix the integral part of its vertex (which fits, the columns being
    nonnegative); the residual shrinks, so dominance and peeling remove candidates, and the search is
    tried again. The CBC integer program is the last resort, when a vertex has no integral part.
    """
    x = {}
    resid = dict(target)
    alive = list(range(len(vecs)))
    while True:
        peeled = _peel(vecs, resid, alive)
        if peeled is None:
            return None
        forced, alive = peeled
        for i, n in forced.items():
            x[i] = x.get(i, 0) + n
        if not resid:
            return x
        try:
            sol = _cover_search(vecs, resid, alive)
        except _SearchBudget:
            sol = None
        else:
            if sol is None:
                return None
            for i, n in sol.items():
                x[i] = x.get(i, 0) + n
            return x
        relaxed = _solve_with_cbc([vecs[i] for i in alive], resid, msg, integer=False)
        if relaxed is None:
            return None
        fixed = {i: int(v + 1e-6) for i, v in zip(alive, relaxed) if v >= 1 - 1e-6}
        if not fixed:
            sol = _solve_with_cbc([vecs[i] for i in alive], resid, msg, integer=True)
            if sol is None:
                return None
            for i, n in zip(alive, sol):
                if n:
                    x[i] = x.get(i, 0) + n
            return x
        for i, n in fixed.items():
            for m, c in vecs[i].items():
                resid[m] -= n * c
            x[i] = x.get(i, 0) + n
        if any(c < 0 for c in resid.values()):
            return None


def _decompose(target, n1, n3, msg):
    """``[(n, pairs)]`` with ``sum n * prod (y_k + z_j) == target`` (homogeneous, nonnegative), or raise."""
    cands = _enumerate_candidates(target, n1, n3)
    x = _solve_positive_system([vec for _, vec in cands], target, msg)
    if x is None:
        raise Exception
    return [(n, cands[i][0]) for i, n in x.items() if n]


class _SearchBudget(Exception):
    pass


_SEARCH_NODE_LIMIT = 5000


def _cover_search(vecs, resid, alive, node_limit=_SEARCH_NODE_LIMIT):
    """Nonnegative integers ``{i: n}`` with ``sum n * vecs[i] == resid`` over the candidates ``alive``, or
    ``None``; raises ``_SearchBudget`` after ``node_limit`` nodes.

    At each node the forced candidates are peeled off, then the monomial of the residual covered by the
    fewest fitting candidates is chosen and one copy of one of them is subtracted: candidate ``t`` is
    used, the candidates listed before it are not (removed), so the branches partition the solutions
    and every branch strictly shrinks the residual. ``resid`` and ``alive`` are not modified.
    """
    nodes = 0

    def dfs(resid, alive):
        nonlocal nodes
        nodes += 1
        if nodes > node_limit:
            raise _SearchBudget
        peeled = _peel(vecs, resid, alive)
        if peeled is None:
            return None
        forced, alive = peeled
        if not resid:
            return forced
        cover = {}
        for i in alive:
            for m in vecs[i]:
                cover.setdefault(m, []).append(i)
        m = min(resid, key=lambda mm: len(cover.get(mm, ())))
        idx = cover.get(m, [])
        for pos, t in enumerate(idx):
            sub = dict(resid)
            for mm, c in vecs[t].items():
                sub[mm] -= c
            rest = [i for i in alive if i not in idx[:pos]]
            found = dfs(sub, rest)
            if found is not None:
                found[t] = found.get(t, 0) + 1
                for i, n in forced.items():
                    found[i] = found.get(i, 0) + n
                return found
        return None

    return dfs(dict(resid), list(alive))


def compute_positive_rep(val, var2=None, var3=None, msg=False):
    """Express ``val`` as a nonnegative-integer combination of product-of-differences monomials.

    ``val`` must be a polynomial in ``var2``/``var3`` known (by positivity of
    double Schubert structure constants) to admit an expansion
    ``sum_b n_b * prod (var2_i - var3_j)`` with ``n_b >= 0`` integers. Builds a
    candidate spanning set of such product monomials from ``val``'s own
    monomials, then solves an integer program (via PuLP) for nonnegative
    integer coefficients ``n_b`` matching ``val`` exactly.

    All polynomial arithmetic is on sparse exponent-tuple dicts with the sign ``(-1)^(z-degree)`` of
    every monomial stripped, which makes every candidate vector and (for a positive ``val``) the target
    nonnegative. Candidates are enumerated with pruning against the target
    (:func:`_enumerate_candidates`), the forced ones are peeled off (:func:`_peel`), and what remains
    is solved by an exact search in Python (:func:`_cover_search`) or, if that exceeds its budget, by
    an LP dive and CBC (:func:`_solve_positive_system`). The result is verified against ``val``'s full
    coefficient vector.

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
    try:
        return int(expand(val))
    except Exception:
        pass
    frees = val.free_symbols
    varsimp2 = sorted([m for m in frees if var2.index(m) != -1], key=var2.index)
    varsimp3 = sorted([m for m in frees if var3.index(m) != -1], key=var3.index)
    n1, n3 = len(varsimp2), len(varsimp3)

    terms = _monomial_dict(val, varsimp2 + varsimp3)
    target = {m: c if sum(m[n1:]) % 2 == 0 else -c for m, c in terms.items()}
    if any(c < 0 for c in target.values()):
        raise Exception  # not a positive combination of products of differences

    # products of differences are homogeneous: each degree is decomposed on its own
    by_degree = {}
    for m, c in target.items():
        by_degree.setdefault(sum(m), {})[m] = c
    solution = []
    for component in by_degree.values():
        solution += _decompose(component, n1, n3, msg)

    total = {}
    for n, pairs in solution:
        vec = {(0,) * (n1 + n3): n}
        for k, j in pairs:
            vec = _vec_mul(vec, _difference_vec((k,), j, n1, n3))
        for m, c in vec.items():
            total[m] = total.get(m, 0) + c
    if total != target:
        raise Exception
    val2 = S.Zero
    for n, pairs in solution:
        val2 += n * prod([varsimp2[k] - varsimp3[j] for k, j in pairs], start=S.One)
    return val2


def _solve_with_cbc(cols, rhs, msg, integer=True):
    """Nonnegative ``x`` with ``sum_i x_i cols[i] == rhs`` by PuLP/CBC: integers, or with ``integer=False``
    a vertex of the LP relaxation (floats). ``None`` if infeasible."""
    import pulp as pu

    lp_prob = pu.LpProblem("Problem", pu.LpMinimize)
    vrs = [lp_prob.add_variable(f"a{i}", lowBound=0, cat="Integer" if integer else "Continuous") for i in range(len(cols))]
    # feasibility is all that is wanted of the integer program (a nonzero objective would make CBC prove
    # optimality); for the relaxation a sparse vertex is preferable
    lp_prob += 0 if integer else pu.lpSum(vrs)
    eqs = {}
    for var, col in zip(vrs, cols):
        for m, c in col.items():
            eqs[m] = eqs.get(m, 0) + (var if c == 1 else c * var)
    for m, lhs in eqs.items():
        lp_prob += lhs == rhs.get(m, 0)
    try:
        solver = cbc_solver(msg)
        status = lp_prob.solve(solver)
    except KeyboardInterrupt:
        current_process = psutil.Process()
        children = current_process.children(recursive=True)
        for child in children:
            child_process = psutil.Process(child.pid)
            child_process.terminate()
            child_process.kill()
        raise KeyboardInterrupt()
    if status != pu.LpStatusOptimal:
        return None
    values = [var.value() or 0.0 for var in vrs]
    if not integer:
        return values
    # solvers return near-integers like 0.9999999999996: round, don't truncate
    return [round(v) for v in values]


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
    # logger.debug(f"NEW {val=} {u2=} {v2=} {w2=}")
    oldval = val
    if u2.inv + v2.inv - w2.inv == 0:
        # logger.debug(f"Hmm this is probably not or val inty true {val=}")
        return val

    if not sign_only and expand(val) == 0:
        # logger.debug(f"Hmm this is probably not true {u2=} {v2=} {w2=} {val=}")
        return 0
    # logger.debug("proceeding")
    u, v, w = u2, v2, w2
    # u, v, w = try_reduce_v(u2, v2, w2)
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
        # logger.debug(f"Return 0 ")
        return 0
    # logger.debug(f"Reduced to {u2=} {v2=} {w2=} {val=}")
    if is_coeff_irreducible(u, v, w):
        u3, v3, w3 = try_reduce_v(u, v, w)
        if not is_coeff_irreducible(u3, v3, w3):
            u, v, w = u3, v3, w3
        else:
            u3, v3, w3 = try_reduce_u(u, v, w)
            if not is_coeff_irreducible(u3, v3, w3):
                u, v, w = u3, v3, w3
    _split_two_b, _split_two = is_split_two(u, v, w)
    # logger.debug("Recording line number")
    if len([i for i in v.code if i != 0]) == 1:
        # logger.debug("Recording line number")
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
        # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        return elem_sym_poly(
            p - r,
            k + p - 1,
            [-var3[i] for i in range(1, n)],
            [-var2[i] for i in hvarset],
        )
        # if expand(val - oldval) != 0:
        #     # logger.debug("This is bad")
        #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
    ## TODO: DOOO THIS
    # cd = code(v)

    # # find one that is hanging off
    # # max_cd = 0
    # # max_index = -1
    # good_index = None
    # for i in range(len(cd)):
    #     if i < len(cd) and cd[i] < cd[i+1]:
    #         continue
    #     good = True
    #     for j in range(i):
    #         if cd[i] <= cd[j] + i - j - 1:
    #             good = False
    #             break
    #     if good:
    #         good_index = i
    #         p = 1
    #         for plus in range(1,len(cd)-i):
    #             if cd[i+plus] < cd[i + plus - 1]:
    #                 break
    #             elif cd[i + plus] > cd[i + plus - 1]:
    #                 break

    #         break
    # if good_index:
    #     if good_index == len(cd) - 1:
    #         p = cd[good_index] - cd[good_index - 1]
    #     else:
    #         p = cd[good_index] - max(cd[good_index+1],cd[good_index-1])
    #     print(f"Broinkspat {u.code=} {v.code=} {w.code=} {good_index=} {p=}")
    #     cd[good_index] -= p
    #     from schubmult.abc import x
    #     from schubmult.rings import DoubleSchubertRing
    #     #from schubmult.symmetric_polynomials import H
    #     R = DoubleSchubertRing(x,var2)
    #     R2 = DoubleSchubertRing(x,var3)
    #     new_v = uncode(cd)
    #     print(f"{code(new_v)=}")
    #     elem = R(u)*R2(new_v)
    #     val2 = S.Zero
    #     R3 = DoubleSchubertRing(x,var3[cd[good_index]+1 - p:])
    #     elem_sym = R3(uncode(([0]*good_index)+[p]))
    #     for w3, vv2 in elem.items():
    #         elem2 = R(w3)*elem_sym
    #         if w in elem2:
    #             val2 += elem2[w] * posify(vv2,u,new_v,w3,var2,var3,msg,sign_only,optimize)
    #     return val2

    # elem_sym_poly(1, good_index + 1, , var3[cd[good_index]+1:])*posify(schubmult_double_pair(u, uncode(cd), var2, var3).get(w,S.Zero),u,uncode(cd),w,var2,var3,msg,sign_only,optimize)

    if will_formula_work(v, u) or u.dominates(w):
        # logger.debug("Recording line number")
        if sign_only:
            return 0
        return dualcoeff(u, v, w, var2, var3)
        # if expand(val - oldval) != 0:
        # logger.debug("This is bad")
        # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=} {will_formula_work(v,u)=} {dominates(u,w)=}")
        # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
    if not v.has_pattern([1, 4, 2, 3]) and not v.has_pattern([4, 1, 3, 2]) and not v.has_pattern([3, 1, 4, 2]) and not v.has_pattern([1, 4, 3, 2]):
        logger.debug("Recording new characterization was used")
        return schubmult_double({u: 1}, v, var2, var3).get(w, 0)

    if w.inv - u.inv == 1:
        # logger.debug("Recording line number")
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
        # if expand(val - oldval) != 0:
        # logger.debug("This is bad")
        # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
        # logger.debug(f"good to go {u2=} {v2=} {w2=}")
        return val
    # if split_two_b:
    #     # logger.debug("Recording line number")
    #     if sign_only:
    #         return 0
    #     cycles = split_two
    #     a1, b1 = cycles[0]
    #     a2, b2 = cycles[1]
    #     a1 -= 1
    #     b1 -= 1
    #     a2 -= 1
    #     b2 -= 1
    #     spo = sorted([a1, b1, a2, b2])
    #     real_a1 = min(spo.index(a1), spo.index(b1))
    #     real_a2 = min(spo.index(a2), spo.index(b2))
    #     real_b1 = max(spo.index(a1), spo.index(b1))
    #     real_b2 = max(spo.index(a2), spo.index(b2))

    #     good1 = False
    #     good2 = False
    #     if real_b1 - real_a1 == 1:
    #         good1 = True
    #     if real_b2 - real_a2 == 1:
    #         good2 = True
    #     a, b = -1, -1
    #     if good1 and not good2:
    #         a, b = min(a2, b2), max(a2, b2)
    #     if good2 and not good1:
    #         a, b = min(a1, b1), max(a1, b1)
    #     arr = [[[], v]]
    #     d = -1
    #     for i in range(len(v) - 1):
    #         if v[i] > v[i + 1]:
    #             d = i + 1
    #     for i in range(d):
    #         arr2 = []

    #         if i in [a1, b1, a2, b2]:
    #             continue
    #         i2 = 1
    #         i2 += len([aa for aa in [a1, b1, a2, b2] if i > aa])
    #         for vr, v2 in arr:
    #             dpret = pull_out_var(i2, v2)
    #             for v3r, v3 in dpret:
    #                 arr2 += [[[*vr, (v3r, i + 1)], v3]]
    #         arr = arr2
    #     val = 0

    #     if good1:
    #         arr2 = []
    #         for L in arr:
    #             v3 = L[-1]
    #             if v3[real_a1] < v3[real_b1]:
    #                 continue
    #             v3 = v3.swap(real_a1, real_b1)
    #             arr2 += [[L[0], v3]]
    #         arr = arr2
    #         if not good2:
    #             for i in range(4):
    #                 arr2 = []

    #                 if i in [real_a2, real_b2]:
    #                     continue
    #                 if i == real_a1:
    #                     var_index = min(a1, b1) + 1
    #                 elif i == real_b1:
    #                     var_index = max(a1, b1) + 1
    #                 i2 = 1
    #                 i2 += len([aa for aa in [real_a2, real_b2] if i > aa])
    #                 for vr, v2 in arr:
    #                     dpret = pull_out_var(i2, v2)
    #                     for v3r, v3 in dpret:
    #                         arr2 += [[[*vr, (v3r, var_index)], v3]]
    #                 arr = arr2
    #     if good2:
    #         arr2 = []
    #         for L in arr:
    #             v3 = L[-1]
    #             try:
    #                 if v3[real_a2] < v3[real_b2]:
    #                     continue
    #                 v3 = v3.swap(real_a2, real_b2)
    #             except IndexError:
    #                 continue
    #             arr2 += [[L[0], v3]]
    #         arr = arr2
    #         if not good1:
    #             for i in range(4):
    #                 arr2 = []

    #                 if i in [real_a1, real_b1]:
    #                     continue
    #                 i2 = 1
    #                 i2 += len([aa for aa in [real_a1, real_b1] if i > aa])
    #                 if i == real_a2:
    #                     var_index = min(a2, b2) + 1
    #                 elif i == real_b2:
    #                     var_index = max(a2, b2) + 1
    #                 for vr, v2 in arr:
    #                     dpret = pull_out_var(i2, v2)
    #                     for v3r, v3 in dpret:
    #                         arr2 += [[[*vr, (v3r, var_index)], v3]]
    #                 arr = arr2

    #         for L in arr:
    #             v3 = L[-1]
    #             tomul = 1
    #             doschubpoly = True
    #             if (not good1 or not good2) and v3[0] < v3[1] and (good1 or good2):
    #                 continue
    #             if (good1 or good2) and (not good1 or not good2):
    #                 v3 = v3.swap(0, 1)
    #             elif not good1 and not good2:
    #                 doschubpoly = False
    #                 if v3[0] < v3[1]:
    #                     dual_u = uncode([2, 0])
    #                     dual_w = Permutation([4, 2, 1, 3])
    #                     coeff = perm_act(dualcoeff(dual_u, v3, dual_w, var2, var3), 2, var2)

    #                 elif len(v3) < 3 or v3[1] < v3[2]:
    #                     if len(v3) <= 3 or v3[2] < v3[3]:
    #                         coeff = 0
    #                         continue
    #                     v3 = v3.swap(0, 1).swap(2, 3)
    #                     coeff = perm_act(schubpoly(v3, var2, var3), 2, var2)
    #                 elif len(v3) <= 3 or v3[2] < v3[3]:
    #                     if len(v3) <= 3:
    #                         v3 += [4]
    #                     v3 = v3.swap(2, 3)
    #                     coeff = perm_act(
    #                         posify(
    #                             schubmult_one(Permutation([1, 3, 2]), v3, var2, var3).get(
    #                                 Permutation([2, 4, 3, 1]),
    #                                 0,
    #                             ),
    #                             Permutation([1, 3, 2]),
    #                             v3,
    #                             Permutation([2, 4, 3, 1]),
    #                             var2,
    #                             var3,
    #                             msg,
    #                             do_pos_neg,
    #                             optimize=optimize,
    #                         ),
    #                         2,
    #                         var2,
    #                     )
    #                     # logger.debug(f"{coeff=}")
    #                 else:
    #                     coeff = perm_act(
    #                         schubmult_one(Permutation([1, 3, 2]), v3, var2, var3).get(
    #                             Permutation([2, 4, 1, 3]),
    #                             0,
    #                         ),
    #                         2,
    #                         var2,
    #                     )
    #                 # logger.debug(f"{coeff=}")
    #                 # if expand(coeff) == 0:
    #                 #     # logger.debug("coeff 0 oh no")
    #                 tomul = sympify(coeff)
    #             toadd = 1
    #             for i in range(len(L[0])):
    #                 var_index = L[0][i][1]
    #                 oaf = L[0][i][0]
    #                 if var_index - 1 >= len(w):
    #                     yv = var_index
    #                 else:
    #                     yv = w[var_index - 1]
    #                 for j in range(len(oaf)):
    #                     toadd *= var2[yv] - var3[oaf[j]]
    #             if (not good1 or not good2) and (good1 or good2):
    #                 varo = [0, var2[w[a]], var2[w[b]]]
    #             else:
    #                 varo = [0, *[var2[w[spo[k]]] for k in range(4)]]
    #             if doschubpoly:
    #                 toadd *= schubpoly(v3, varo, var3)
    #             else:
    #                 subs_dict3 = {var2[i]: varo[i] for i in range(len(varo))}
    #                 toadd *= efficient_subs(tomul, subs_dict3)
    #             val += toadd
    #             # logger.debug(f"accum {val=}")
    #         #logger.debug(f"{expand(val-oldval)=}")
    #         # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
    #         return val
    if will_formula_work(u, v):
        # logger.debug("Recording line number")
        if sign_only:
            return 0
        # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        return forwardcoeff(u, v, w, var2, var3)
        # if expand(val - oldval) != 0:
        #     # logger.debug("This is bad")
        #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
    # logger.debug("Recording line number")
    # c01 = code(u)
    # c02 = code(w)
    # c03 = code(v)

    c1 = (~u).code
    c2 = (~w).code

    if u.one_dominates(w):
        if sign_only:
            return 0
        while c1[0] != c2[0]:
            w = w.swap(c2[0] - 1, c2[0])
            v = v.swap(c2[0] - 1, c2[0])
            # w[c2[0] - 1], w[c2[0]] = w[c2[0]], w[c2[0] - 1]
            # v[c2[0] - 1], v[c2[0]] = v[c2[0]], v[c2[0] - 1]
            # w = tuple(w)
            # v = tuple(v)
            c2 = (~w).code
            # c03 = code(v)
            # c01 = code(u)
            # c02 = code(w)
        # if is_reducible(v):
        #     # logger.debug("Recording line number")
        #     if sign_only:
        #         return 0
        #     newc = []
        #     elemc = []
        #     for i in range(len(c03)):
        #         if c03[i] > 0:
        #             newc += [c03[i] - 1]
        #             elemc += [1]
        #         else:
        #             break
        #     v3 = uncode(newc)
        #     coeff_dict = schubmult_one(
        #         u,
        #         uncode(elemc),
        #         var2,
        #         var3,
        #     )
        #     val = 0
        #     for new_w in coeff_dict:
        #         tomul = coeff_dict[new_w]
        #         newval = schubmult_one(new_w, uncode(newc), var2, var3).get(
        #             w,
        #             0,
        #         )
        #         # logger.debug(f"Calling posify on {newval=} {new_w=} {uncode(newc)=} {w=}")
        #         newval = posify(newval, new_w, uncode(newc), w, var2, var3, msg, do_pos_neg, optimize=optimize,elem_sym=elem_sym)
        #         val += tomul * shiftsubz(newval)
        #     # if expand(val - oldval) != 0:
        #     #     # logger.debug("This is bad")
        #     #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
        #         # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        #     return val
        # removed, iffy (hard to implement)
        # if c01[0] == c02[0] and c01[0] != 0:
        #     # logger.debug("Recording line number")
        #     if sign_only:
        #         return 0
        #     varl = c01[0]
        #     u3 = uncode([0] + c01[1:])
        #     w3 = uncode([0] + c02[1:])
        #     val = 0
        #     val = schubmult_one(u3, v, var2, var3).get(
        #         w3,
        #         0,
        #     )
        #     # logger.debug(f"Calling posify on {val=} {u3=} {v=} {w3=}")
        #     val = posify(val, u3, v, w3, var2, var3, msg, do_pos_neg, optimize=optimize,elem_sym=elem_sym)
        #     for i in range(varl):
        #         val = perm_act(val, i + 1, var2)
        #     # if expand(val - oldval) != 0:
        #     #     # logger.debug("This is bad")
        #     #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
        #     # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        #     return val
        if c1[0] == c2[0]:
            # logger.debug("Recording line number")
            if sign_only:
                return 0
            vp = pull_out_var(c1[0] + 1, v)
            u3 = phi1(u)
            w3 = phi1(w)
            val = 0
            for arr, v3 in vp:
                tomul = 1
                for i in range(len(arr)):
                    tomul *= var2[1] - var3[arr[i]]

                val2 = schubmult_double_pair(u3, v3, var2, var3).get(
                    w3,
                    0,
                )
                val2 = posify(val2, u3, v3, w3, var2, var3, msg, optimize=optimize)
                val += tomul * shiftsub(val2, var2)
            # if expand(val - oldval) != 0:
            #     # logger.debug("This is bad")
            #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=")
            # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
            return val
    # logger.debug("Fell all the way through. Cleverness did not save us")
    if not sign_only:
        # logger.debug("Recording line number")
        if optimize:
            # if elem_sym:
            #     # print(f"{elem_sym=}")
            #     val2 = compute_positive_rep_new(elem_sym, var2, var3, msg, False)
            if u.inv + v.inv - w.inv == 1:
                val2 = compute_positive_rep(val, var2, var3, msg)
            else:
                val2 = compute_positive_rep(val, var2, var3, msg)
            if val2 is not None:
                val = val2
            return val
        if optimize is None:
            raise Exception
        # logger.debug("RETURNINGOLDVAL")
        return oldval
    # logger.debug("Recording line number")
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
        # logger.debug("Warning, failed on a case")
        raise Exception(f"{val=} {val2=} {u2=} {v2=} {w2=}")
    # print("FROFL")
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
        # logger.debug(f"{coeff_dict.get(w,0)=} {w=} {perm=} {vmun1=} {v=} {muv=}")
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
        # logger.debug("Recording line number")
        vp = v * (~perm)
        if vp.inv == v.inv - perm.inv:
            return schubpoly(vp, var2, var3)
    dpret = []
    ret = 0
    if u.dominates(perm):
        dpret = dualpieri(u, v, perm)
    else:
        # logger.debug("Recording line number")
        dpret = []
        # logger.debug("Recording line number")
        th = u.theta()
        muu = uncode(th)
        umun1 = (~u) * muu
        w = perm * umun1
        # logger.debug("spiggle")
        # logger.debug(f"{u=} {muu=} {v=} {w=} {perm=}")
        # logger.debug(f"{w=} {perm=}")
        if w.inv == umun1.inv + perm.inv:
            dpret = dualpieri(muu, v, w)
            # logger.debug(f"{muu=} {v=} {w=}")
            # logger.debug(f"{dpret=}")
    for vlist, vp in dpret:
        # logger.debug("Recording line number")
        toadd = 1
        for i in range(len(vlist)):
            for j in range(len(vlist[i])):
                toadd *= var2[i + 1] - var3[vlist[i][j]]
        toadd *= schubpoly(vp, var2, var3, len(vlist) + 1)
        ret += toadd
    return ret
    # logger.debug("Recording line number")
    # schub_val = schubmult_one(u, v, var2, var3)
    # val_ret = schub_val.get(perm, 0)
    # if expand(val - val_ret) != 0:
    #     # logger.debug(f"{schub_val=}")
    #     # logger.debug(f"{val=} {u=} {v=} {var2[1]=} {var3[1]=}  {perm=} {schub_val.get(perm,0)=}")
    # logger.debug(f"good to go {ret=}")


def dualpieri(mu, v, w):
    """Dual Pieri expansion used by ``dualcoeff``: enumerate the data witnessing
    ``S_mu * S_v -> S_w`` when ``mu`` is dominant.

    Compares ``mu``'s inverse code against ``w``'s inverse code layer by layer,
    peeling one "cycle" of variables per layer via ``divdiffable``/``pull_out_var``,
    and returns the list of ``[vlist, vp]`` pairs consumed by ``dualcoeff`` to
    build the final positive expression (empty list if ``w`` is not reachable
    from ``mu``, ``v`` this way).
    """
    # logger.debug(f"dualpieri {mu=} {v=} {w=}")
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
    # c = permtrim(c)
    # logger.debug("Recording line number")
    res = [[[], v]]
    # logger.debug(f"{v=} {type(v)=}")
    for i in range(len(lm)):
        # logger.debug(f"{res=}")
        res2 = []
        for vlist, vplist in res:
            vp = vplist
            vpl = divdiffable(vp, cycle(lm[i] + 1, cn1w[i] - lm[i]))
            # logger.debug(f"{vpl=} {type(vpl)=}")
            if len(vpl) == 0:
                continue
            vl = pull_out_var(lm[i] + 1, vpl)
            # logger.debug(f"{vl=}")
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
    # logger.debug(f"{res2=}")
    return res2
