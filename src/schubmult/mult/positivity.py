"""Manifestly positive representations of double Schubert and Grothendieck structure constants.

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
"""

from __future__ import annotations

import concurrent.futures as cf
import operator
import shutil
import tempfile
import threading
from collections import deque
from functools import cache
from typing import TYPE_CHECKING

import psutil
from cachetools import cached
from cachetools.keys import hashkey

from schubmult.abc import beta
from schubmult.combinatorics.permutation import (
    Permutation,
    cycle,
    phi1,
    uncode,
)
from schubmult.symbolic import Add, Mul, Pow, S, Symbol, expand, prod, sympify, sympify_sympy
from schubmult.symbolic.common_polys import _vars, efficient_subs, elem_sym_poly, schubpoly
from schubmult.utils.logging import get_logger
from schubmult.utils.schub_lib import (
    divdiffable,
    is_coeff_irreducible,
    pull_out_var,
    reduce_coeff,
    reduce_descents,
    try_reduce_u,
    try_reduce_v,
    will_formula_work,
)

from .double import schubmult_double, schubmult_double_pair, schubmult_double_pair_generic_alt

if TYPE_CHECKING:
    from schubmult._typing import Alphabet

logger = get_logger(__name__)

_posify_cache = {}

# graceful failure ceiling: candidate enumeration and the LP grow combinatorially,
# and unbounded growth has OOM-killed the whole machine before
_MAX_CANDIDATES = 1_000_000


def _var_index(varlist: Alphabet, sym: Symbol) -> int:
    """Position of ``sym`` in ``varlist`` or ``-1``.

    ``posify`` recurses with a shifted slice of ``var2``, turning a ``GeneratingSet`` into a plain tuple
    whose ``.index`` raises instead of returning ``-1`` for foreign symbols.
    """
    try:
        return varlist.index(sym)
    except ValueError:
        return -1


def cbc_solver(msg: bool = False):
    """PuLP ``COIN_CMD`` using the CBC binary bundled by ``pulp[cbc]``, falling back to ``cbc`` on PATH."""
    import pulp as pu

    try:
        from cbcbox import cbc_bin_path

        return pu.COIN_CMD(msg=msg, path=cbc_bin_path())
    except ImportError:
        return pu.COIN_CMD(msg=msg)


def _solver(msg, tmpdir=None):
    """The MILP solver: HiGHS if available, else CBC.  ``tmpdir`` receives its temporary files."""
    import pulp as pu

    if pu.HiGHS_CMD().available():
        # HiGHS's presolve costs more than it saves on these equality systems (they are already reduced by _presolve).
        # One thread per solver: concurrency comes from running several solvers, and the default thread pool of
        # HiGHS only spins, which costs about 40% more CPU time without making a solve any faster
        solver = pu.HiGHS_CMD(msg=msg, options=["presolve=off", "threads=1"])
    else:
        solver = cbc_solver(msg)
    if tmpdir is not None:
        solver.tmpDir = tmpdir
    return solver


def _is_optimal(status):
    """``LpProblem.solve`` returns an int status in PuLP < 4 and an ``LpSolveStats`` object in PuLP >= 4."""
    if status is None:
        return False
    return int(getattr(status, "status", status)) == 1


def _exact_int(x):
    """``x`` as a Python ``int`` if it is an exact integer (Python or symbolic), else ``None``."""
    if isinstance(x, int):
        return x
    return int(x) if getattr(x, "is_Integer", False) else None


def _weighted_degree(expr, weights):
    """Degree of ``expr`` under the grading ``weights`` (``{symbol: weight}``), or ``None`` if not provably homogeneous.

    Purely structural, so it never expands: every sum must have summands of one common degree.
    """
    if expr.is_Number:
        return 0
    if expr.is_Symbol:
        return weights.get(expr)
    if expr.is_Add:
        degs = {_weighted_degree(arg, weights) for arg in expr.args}
        return degs.pop() if len(degs) == 1 and None not in degs else None
    if expr.is_Mul:
        total = 0
        for arg in expr.args:
            deg = _weighted_degree(arg, weights)
            if deg is None:
                return None
            total += deg
        return total
    if expr.is_Pow and expr.args[1].is_Integer:
        deg = _weighted_degree(expr.args[0], weights)
        return None if deg is None else deg * int(expr.args[1])
    return None


def _raise_if_cancelled(stop):
    """Raise ``RuntimeError`` if the ``stop`` event of an abandoned run (or ``None``) is set; worker threads cannot be killed, so they poll."""
    if stop is not None and stop.is_set():
        raise RuntimeError("integer program cancelled")


def _presolve(cands, target, stop=None):
    """Exact presolve of ``sum_j a_j * cands[j] == target`` over integers ``a_j >= 0``.

    ``cands`` is a list of ``{row: coefficient}`` dicts.  Repeatedly applies rules that are implied by
    non-negativity and integrality: singleton rows fix their column, zero rows with one-signed columns
    zero them, and a row with a positive (negative) target and no negative (positive) coefficient bounds
    every column in it by the target.  Returns ``None`` if the system is infeasible, otherwise
    ``(fixed, active, resid)``: the values fixed so far, the columns still free and the residual target.
    A ``stop`` event that gets set cancels the presolve with a ``RuntimeError``.
    """
    rows = {}
    for j, svec in enumerate(cands):
        _raise_if_cancelled(stop)
        for k in svec:
            rows.setdefault(k, set()).add(j)
    resid = dict(target)
    active = set(range(len(cands)))
    fixed = {}
    queue = deque(set(rows) | set(resid))
    queued = set(queue)

    def remove(j, x):
        active.discard(j)
        if x:
            fixed[j] = x
        for k2, c2 in cands[j].items():
            if x:
                resid[k2] = resid.get(k2, 0) - x * c2
            rows[k2].discard(j)
            if k2 not in queued:
                queued.add(k2)
                queue.append(k2)

    while queue:
        _raise_if_cancelled(stop)
        k = queue.popleft()
        queued.discard(k)
        t = resid.get(k, 0)
        cols = rows.get(k)
        if not cols:
            if t != 0:
                return None
            continue
        if len(cols) == 1:
            (j,) = cols
            c = cands[j][k]
            if t % c or t // c < 0:
                return None
            remove(j, t // c)
            continue
        pos = [j for j in cols if cands[j][k] > 0]
        neg = len(cols) - len(pos)
        if (t > 0 and not pos) or (t < 0 and not neg):
            return None
        if t == 0 and (not pos or not neg):
            for j in list(cols):
                remove(j, 0)
        elif t > 0 and not neg:
            for j in [j for j in cols if cands[j][k] > t]:
                remove(j, 0)
        elif t < 0 and not pos:
            for j in [j for j in cols if -cands[j][k] > -t]:
                remove(j, 0)
    return fixed, active, resid


def _solve_milp(svecs, target, msg=False, tmpdir=None, stop=None):
    """Non-negative integer multiplicities ``a`` with ``sum_j a[j] * svecs[j] == target``, or ``None`` if there are none.

    ``svecs`` is a list of ``{row: coefficient}`` dicts.  Only integer arithmetic and a solver subprocess are involved,
    so concurrent calls are safe.  The solver's temporary files go to ``tmpdir``.  A ``stop`` event that gets set cancels
    the solve with a ``RuntimeError`` as soon as the next check is reached: during the presolve, the model construction,
    and just before the solver is launched.
    """
    import pulp as pu

    pre = _presolve(svecs, target, stop)
    if pre is None:
        return None
    fixed, active, resid = pre
    values = [0] * len(svecs)
    for j, x in fixed.items():
        values[j] = x
    if not active:
        return values
    cols = sorted(active)
    lp_prob = pu.LpProblem("Problem", pu.LpMinimize)
    vrs = {}
    for j in cols:
        _raise_if_cancelled(stop)
        vrs[j] = lp_prob.add_variable(f"a{j}", lowBound=0, cat="Integer")
    # minimizing the total multiplicity solves about twice as fast as a feasibility-only objective
    lp_prob += pu.lpSum(vrs.values())
    eqs = {}
    for j in cols:
        _raise_if_cancelled(stop)
        for k, c in svecs[j].items():
            eqs.setdefault(k, {})[vrs[j]] = c
    for k, terms in eqs.items():
        _raise_if_cancelled(stop)
        # lpSum rather than LpAffineExpression(dict): the latter cannot form constraints in PuLP >= 4
        lp_prob += pu.lpSum(c * v for v, c in terms.items()) == resid.get(k, 0)
    _raise_if_cancelled(stop)
    status = lp_prob.solve(_solver(msg, tmpdir))
    if not _is_optimal(status):
        return None
    # read values here: in PuLP >= 4 variables become unusable once lp_prob is garbage collected
    for j, v in vrs.items():
        x = v.value()
        # round, don't truncate: solvers return near-integers like 0.9999999999996
        values[j] = round(x) if x is not None else 0
    return values


def _resolve_workers(max_workers, n_items):
    """Number of concurrent solves for ``n_items`` coefficients.

    ``max_workers`` if given, otherwise one.
    """
    if max_workers is None:
        workers = 1
    else:
        workers = operator.index(max_workers)
        if workers < 1:
            raise ValueError(f"max_workers must be at least 1, got {max_workers}")
    return max(1, min(workers, n_items))


def _terminate_solvers(tmpdir):
    """Kill the solver processes working in ``tmpdir``, which show up as the temporary files on their command lines."""
    for child in psutil.Process().children(recursive=True):
        try:
            if any(tmpdir in arg for arg in child.cmdline()):
                child.kill()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue


def _abort_solves(futures, tmpdir, stop):
    """Abandon the integer programs ``futures`` of an interrupted run: cancel the queued, kill the solvers of the running."""
    stop.set()
    for fut in futures:
        fut.cancel()
    live = [fut for fut in futures if not fut.done()]
    while True:
        # Workers that are still building a model notice ``stop`` within moments, so give them the chance to finish before
        # scanning the process table: every file of the scan gives up the GIL, and next to threads that are busy in Python
        # it takes a full switch interval to get it back, which makes the scan crawl for seconds.  What remains is waiting
        # for a solver.  The loop repeats, since a worker that was building its model during an earlier sweep launches its
        # solver later.
        live = list(cf.wait(live, timeout=0.1).not_done)
        # at least one sweep, for the solver of an inline solve that was interrupted before it had a future
        _terminate_solvers(tmpdir)
        if not live:
            return


class _InlineExecutor:
    """Runs each job at once on the calling thread, for a single worker that needs no threads."""

    def submit(self, fn, *args):
        fut = cf.Future()
        try:
            fut.set_result(fn(*args))
        except Exception as exc:
            fut.set_exception(exc)
        return fut

    def shutdown(self, wait=True):
        pass


def _posify_steps(val, var2, var3):
    """Coroutine behind :func:`groth_posify`: all the symbolic work, leaving the integer programs to the caller.

    Yields each integer program as ``(svecs, target)``, to be answered with ``send`` of its :func:`_solve_milp`
    solution (``None`` if there is none).  The result of ``groth_posify`` is the return value.
    """
    import itertools

    try:
        return int(val)
    except Exception:
        pass
    bs = sympify(beta)
    e = sympify(val)

    # homogeneous degree: yz-degree minus beta-degree is constant across terms
    num0 = expand(sympify(sympify_sympy(e).as_numer_denom()[0]))
    t0 = next(mono for mono, c in num0.as_coefficients_dict().items() if c != 0)
    t0d = t0.as_powers_dict()
    d = sum(int(p) for g, p in t0d.items() if g.is_Symbol and sympify(g) != bs) - int(t0d.get(bs, 0))

    e1 = expand(efficient_subs(e, {bs: S.One}))
    frees = e1.free_symbols
    ysyms = sorted([s for s in frees if var2.index(s) != -1], key=lambda s: var2.index(s))
    zsyms = sorted([s for s in frees if var3.index(s) != -1], key=lambda s: var3.index(s))

    uu = {ys: Symbol(f"u_{var2.index(ys)}") for ys in ysyms}
    vv = {zs: Symbol(f"v_{var3.index(zs)}") for zs in zsyms}
    back = {u: S.One + bs * sympify(ys) for ys, u in uu.items()}
    back.update({v: S.One + bs * sympify(zs) for zs, v in vv.items()})

    fwd = {sympify(ys): uu[ys] - S.One for ys in ysyms}
    fwd.update({sympify(zs): vv[zs] - S.One for zs in zsyms})
    laurent = expand(e1.subs(fwd))

    # clear the (monomial) denominator structurally; everything is a polynomial
    # in u, v from here on and terms are opaque basis vectors
    num_l, den_l = laurent.as_numer_denom()
    num_l = expand(num_l)
    caps = {g: int(p) for g, p in den_l.as_powers_dict().items() if g != S.One and int(p) != 0}

    # packed-int monomial representation over the u, v generators: PBITS bits
    # of exponent per generator, so monomial products are plain integer sums
    PBITS = 10
    gens = [*[uu[ys] for ys in ysyms], *[vv[zs] for zs in zsyms]]
    gshift = {g: PBITS * i for i, g in enumerate(gens)}

    def to_vec(expr):
        vec = {}
        degs = {}
        for mono, c in expr.as_coefficients_dict().items():
            if c == 0:
                continue
            ci = _exact_int(c)
            if ci is None:
                raise ValueError(f"non-integer coefficient {c} in {val}")
            key = 0
            deg = 0
            for g, p in mono.as_powers_dict().items():
                if g in gshift:
                    pi = _exact_int(p)
                    if pi is None or pi < 0:
                        raise ValueError(f"unexpected exponent {p} of {g}")
                    key += pi << gshift[g]
                    deg += pi
                elif g != S.One:
                    raise ValueError(f"unexpected generator {g}")
            vec[key] = vec.get(key, 0) + ci
            degs[key] = deg
        return vec, degs

    def dict_mul(a, b):
        out = {}
        for m1, c1 in a.items():
            for m2, c2 in b.items():
                k = m1 + m2
                cc = out.get(k, 0) + c1 * c2
                if cc:
                    out[k] = cc
                elif k in out:
                    del out[k]
        return out

    vec0, deg0 = to_vec(num_l)
    max_m = max(deg0.values(), default=0)
    if max_m >= 1 << PBITS:
        raise ValueError(f"degree {max_m} too large for the packed monomial representation")
    # candidates are homogeneous of total degree #pairs + clearing-shift degree,
    # so the problem splits into independent blocks by degree and degrees
    # absent from the target can be pruned losslessly
    target_degs = frozenset(deg0.values())
    blocks = {}
    for key, c in vec0.items():
        blocks.setdefault(deg0[key], {})[key] = c

    # construction state is shared across coefficients: difference products
    # depend only on the generators, levels also on caps and target degrees;
    # both are built lazily on demand
    gens_key = (tuple(str(s) for s in ysyms), tuple(str(s) for s in zsyms))
    prods_state = _posify_cache.get(gens_key)
    if prods_state is None:
        # difference products built incrementally: each level multiplies by one
        # 2-term pair vector (integer convolution, no symbolic expand)
        pair_vecs = [({1 << gshift[vv[zs]]: 1, 1 << gshift[uu[ys]]: -1}, (ys, zs)) for ys in ysyms for zs in zsyms]
        prods_state = {"pair_vecs": pair_vecs, "prods": [[({0: 1}, (), 0)]]}
        _posify_cache[gens_key] = prods_state

    cache_key = (gens_key, tuple(sorted((str(g), p) for g, p in caps.items())), target_degs)
    state = _posify_cache.get(cache_key)
    if state is None:
        # clearing monomials as precomputed packed shifts with their degrees
        usage_shifts = []
        for usage in itertools.product(*[range(caps[g] + 1) for g in caps]):
            shift = 0
            sdeg = 0
            for g, a in zip(caps, usage):
                shift += (caps[g] - a) << gshift[g]
                sdeg += caps[g] - a
            usage_shifts.append((shift, dict(zip(caps, usage)), sdeg))
        state = {"usage_shifts": usage_shifts, "levels": {}}
        _posify_cache[cache_key] = state

    def get_level(m, degree):
        """Candidates with ``m`` difference factors and total degree ``degree``."""
        level = state["levels"].get(m)
        if level is None:
            prods = prods_state["prods"]
            pair_vecs = prods_state["pair_vecs"]
            shifts = {}
            for shift, usage, sdeg in state["usage_shifts"]:
                if m + sdeg in target_degs:
                    shifts.setdefault(m + sdeg, []).append((shift, usage))
            level = {}
            if shifts:
                while len(prods) <= m:
                    cur = []
                    for vec, combo, start in prods[-1]:
                        for idx in range(start, len(pair_vecs)):
                            pv, pair = pair_vecs[idx]
                            cur.append((dict_mul(vec, pv), (*combo, pair), idx))
                    prods.append(cur)
                for deg, deg_shifts in shifts.items():
                    level[deg] = [({k + shift: c for k, c in vec.items()}, combo, usage, m) for vec, combo, _ in prods[m] for shift, usage in deg_shifts]
            state["levels"][m] = level
        return level.get(degree, [])

    # solve each degree block separately, escalating the factor-count ceiling:
    # small LPs solve fast and usually suffice
    chosen = []
    for degree, target in sorted(blocks.items()):
        candidates = []
        values = None
        for m in range(max(d, 0), degree + 1):
            candidates = [*candidates, *get_level(m, degree)]
            if len(candidates) > _MAX_CANDIDATES:
                raise ValueError(f"candidate set too large ({len(candidates)}) for {val}")
            values = yield ([svec for svec, _, _, _ in candidates], target)
            if values is not None:
                break
        if values is None:
            raise ValueError(f"no positive representation found for {val}")
        chosen.extend((cand, x) for cand, x in zip(candidates, values) if x)

    result = S.Zero
    total = {}
    for (svec, combo, usage, m), xi in chosen:
        term = S.One
        for ys, zs in combo:
            term = term * (sympify(zs) - sympify(ys))
        for g, a in usage.items():
            if a:
                term = term / back[g] ** a
        result += xi * bs ** (m - d) * term
        for k, c in svec.items():
            total[k] = total.get(k, 0) + xi * c

    # exact verification.  The integer identity is the claim at beta = 1; both sides are homogeneous of
    # degree d (deg y = deg z = 1, deg beta = -1, Laurent atoms 1 + beta*y of degree 0), which determines
    # them from their beta = 1 specialization.  When homogeneity of the input cannot be shown
    # structurally, fall back to expanding the difference.
    weights = {bs: -1, **dict.fromkeys(map(sympify, (*ysyms, *zsyms)), 1)}
    if _weighted_degree(e, weights) == d:
        if {k: c for k, c in total.items() if c} != vec0:
            raise ValueError(f"positive representation check failed for {val}: got {result}")
    else:
        diff_num, _ = sympify_sympy(result - e).as_numer_denom()
        if expand(sympify(diff_num)) != S.Zero:
            raise ValueError(f"positive representation check failed for {val}: got {result}")
    return result


def groth_posify(val, var2, var3, msg=False):
    """Positive FGL form of a coefficient, at ``beta = 1``.

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
    """
    (result,) = groth_posify_many([val], var2, var3, msg, max_workers=1)
    return result


def groth_posify_many(vals, var2, var3, msg=False, *, max_workers=None, progress=None):
    """:func:`groth_posify` of each coefficient in ``vals``, solving the integer programs concurrently.

    The solver subprocesses are driven by up to ``max_workers`` threads, while the symbolic work and the construction
    caches shared between coefficients stay on the calling thread.  The results do not depend on ``max_workers``: they
    are those of sequential :func:`groth_posify` calls, yielded in the order of ``vals`` as soon as each is ready.  If a
    coefficient has no positive representation, the results before it are yielded and then its exception is raised.
    Closing the generator early, or interrupting it, stops the solver processes and removes their temporary files.

    Args:
        vals: Coefficients, as for :func:`groth_posify`.
        var2: First secondary alphabet (``y``).
        var3: Second secondary alphabet (``z``).
        msg: Passed through to the solver as its verbosity option.
        max_workers: Number of integer programs solved at once, at least 1; 1 runs everything on the calling
            thread.  Defaults to 1.  If ``msg`` is enabled with multiple workers, solver logs may interleave.
            Each solver process needs memory of its own, which is large for the hardest cases.
        progress: Called as ``progress(i)``, on the calling thread, when work on ``vals[i]`` starts.

    Returns:
        A generator of the ``groth_posify`` results, in the order of ``vals``.

    Raises:
        ValueError: If ``max_workers`` is less than 1.
    """
    vals = list(vals)
    workers = _resolve_workers(max_workers, len(vals))
    return _posify_many(vals, var2, var3, msg, workers, progress)


def _posify_many(vals, var2, var3, msg, workers, progress):
    """Scheduler of :func:`groth_posify_many`: runs one :func:`_posify_steps` coroutine per coefficient.

    A coroutine is resumed on this thread whenever the integer program it waits for is done.  Coroutines are started in
    order, and the results are released in order.
    """
    n = len(vals)
    tmpdir = tempfile.mkdtemp(prefix="schubmult-posify-")
    stop = threading.Event()
    pool = cf.ThreadPoolExecutor(workers, thread_name_prefix="posify") if workers > 1 else _InlineExecutor()
    # more coroutines than workers, so that the pool keeps busy while this thread prepares the next programs
    limit = 2 * workers if workers > 1 else 1
    gens = {}  # index -> coroutine that has not finished
    futs = {}  # index -> its pending integer program
    results = {}
    failure = None  # (index, exception) of the earliest failed coefficient
    next_start = next_out = 0

    def step(i, fut=None):
        """Resume coroutine ``i`` with the outcome ``fut`` of its last program (none on the first step), and record what it does next."""
        nonlocal failure
        try:
            if fut is None:
                svecs, target = next(gens[i])
            elif fut.exception() is None:
                svecs, target = gens[i].send(fut.result())
            else:
                svecs, target = gens[i].throw(fut.exception())
        except StopIteration as done:
            del gens[i]
            results[i] = done.value
        except Exception as exc:
            del gens[i]
            if failure is None or i < failure[0]:
                failure = (i, exc)
        else:
            futs[i] = pool.submit(_solve_milp, svecs, target, msg, tmpdir, stop)

    try:
        while True:
            while next_out in results:
                yield results.pop(next_out)
                next_out += 1
            if failure is not None and failure[0] == next_out:
                raise failure[1]
            if next_out == n:
                return
            # coefficients after the failed one are never needed
            wanted = [i for i in futs if failure is None or i < failure[0]]
            ready = [i for i in wanted if futs[i].done()]
            if ready:
                for i in ready:
                    step(i, futs.pop(i))
            elif next_start < n and len(gens) < limit and failure is None:
                if progress is not None:
                    progress(next_start)
                gens[next_start] = _posify_steps(vals[next_start], var2, var3)
                step(next_start)
                next_start += 1
            else:
                cf.wait([futs[i] for i in wanted], return_when=cf.FIRST_COMPLETED)
    finally:
        if gens:
            _abort_solves(list(futs.values()), tmpdir, stop)
        pool.shutdown(wait=True)
        shutil.rmtree(tmpdir, ignore_errors=True)


def compute_positive_rep(val: Add, var2: Alphabet, var3: Alphabet, msg: bool = False):
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

    varsimp2 = [m for m in frees if _var_index(var2, m) != -1]
    varsimp3 = [m for m in frees if _var_index(var3, m) != -1]
    varsimp2.sort(key=lambda k: _var_index(var2, k))
    varsimp3.sort(key=lambda k: _var_index(var3, k))

    var22 = tuple(varsimp2)
    var33 = tuple(varsimp3)

    n1 = len(varsimp2)

    base_vectors = {}

    val_expr = expand(val)
    vec0 = {k: v for k, v in val_expr.subs({var3[1]: S.Zero}).as_coefficients_dict().items() if v != S.Zero}
    # val_poly = sympy_poly(val_expr, *var22, *var33)
    # mn = val_poly.monoms()

    def build_single_monom(key, varlist):
        if key.is_Number:
            return tuple([0 for _ in varlist])
        if isinstance(key, Symbol):
            return tuple([1 if key == v else 0 for v in varlist])
        if isinstance(key, Mul):
            result = [0 for _ in varlist]
            for arg in key.args:
                single = build_single_monom(arg, varlist)
                result = [r + s for r, s in zip(result, single)]
            return tuple(result)
        if isinstance(key, Pow):
            base, exp = key.args
            single = build_single_monom(base, varlist)
            return tuple([s * exp for s in single])
        raise ValueError(f"Unsupported key type: {type(key)}")

    def build_monoms(arr, varlist):
        return {build_single_monom(a, varlist) for a in arr}

    mn = build_monoms(Add.make_args(val_expr), [*var22, *var33])
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
    key=lambda val, u2, v2, w2, var2=None, var3=None, msg=False, optimize=True, n=_vars.n: hashkey(val, u2, v2, w2, var2, var3, msg, optimize, n),
)
def posify(
    val,
    u2: Permutation,
    v2: Permutation,
    w2: Permutation,
    var2: Alphabet,
    var3: Alphabet,
    msg=False,
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

    Results are cached by all arguments.

    Args:
        val: The structure constant to re-express positively.
        u2: First factor's permutation.
        v2: Second factor's permutation.
        w2: Target permutation (coefficient of ``S_{w2}``).
        var2: First secondary alphabet.
        var3: Second secondary alphabet.
        msg: Verbosity flag passed down to ``compute_positive_rep``'s LP solver.
        optimize: If ``False``, return ``val`` unchanged when no closed-form
            reduction applies (skip the LP fallback); if ``None`` and that
            case is reached, raise.
        n: Size of the ambient alphabet used when no other bound is available.

    Returns:
        A manifestly positive expression equal to ``val``, or ``val`` itself
        when no formula applies and optimization is disabled.
    """
    if not v2.has_pattern([1, 4, 2, 3]) and not v2.has_pattern([4, 1, 3, 2]) and not v2.has_pattern([3, 1, 4, 2]) and not v2.has_pattern([1, 4, 3, 2]):
        logger.debug("Recording new characterization was used")
        return schubmult_double({u2: 1}, v2, var2, var3).get(w2, 0)
    if u2.inv + v2.inv - w2.inv == 0:
        return val

    if expand(val) == 0:
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

    if is_coeff_irreducible(u, v, w):
        u3, v3, w3 = try_reduce_v(u, v, w)
        if not is_coeff_irreducible(u3, v3, w3):
            u, v, w = u3, v3, w3
        else:
            u3, v3, w3 = try_reduce_u(u, v, w)
            if not is_coeff_irreducible(u3, v3, w3):
                u, v, w = u3, v3, w3

    if len([i for i in v.code if i != 0]) == 1:
        cv = v.code
        k: int = 0
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
        return dualcoeff(u, v, w, var2, var3)

    if not v.has_pattern([1, 4, 2, 3]) and not v.has_pattern([4, 1, 3, 2]) and not v.has_pattern([3, 1, 4, 2]) and not v.has_pattern([1, 4, 3, 2]):
        logger.debug("Recording new characterization was used")
        return schubmult_double({u: 1}, v, var2, var3).get(w, 0)

    if w.inv - u.inv == 1:
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
        return forwardcoeff(u, v, w, var2, var3)

    c1 = (~u).code
    c2 = (~w).code

    if u.one_dominates(w):
        while c1[0] != c2[0]:
            w = w.swap(c2[0] - 1, c2[0])
            v = v.swap(c2[0] - 1, c2[0])

            c2 = (~w).code

        vp = pull_out_var(c1[0] + 1, v)
        u3 = phi1(u)
        w3 = phi1(w)
        val = 0
        for arr, v3 in vp:
            tomul = prod([var2[1] - var3[arr[i]] for i in range(len(arr))])

            # Include the full possible product support, not just the cached prefix.
            shifted_var2 = var2[1 : max(len(u3) + len(v3), len(w3)) + 2]
            val2 = schubmult_double_pair(u3, v3, shifted_var2, var3).get(
                w3,
                0,
            )
            val2 = posify(val2, u3, v3, w3, shifted_var2, var3, msg, optimize=optimize, n=n)
            val += tomul * val2

        return val

    if optimize:
        return compute_positive_rep(val, var2, var3, msg)
    if optimize is None:
        raise ValueError("Optimize is None but no positive formula applied")
    return val


def shiftsub(pol, var2):
    """Shift every ``var2[i]`` in ``pol`` up to ``var2[i + 1]`` (for ``i`` in ``0..98``)."""
    subs_dict = {var2[i]: var2[i + 1] for i in range(99)}
    return efficient_subs(sympify(pol), subs_dict)


def posify_generic_partial(val, u2: Permutation, v2: Permutation, w2: Permutation):
    """``posify`` specialized to the fixed generic alphabets ``_vars.var_g1``/``_vars.var_g2``.

    Asserts (raises on mismatch) that the recomputed positive expression equals
    the input ``val``, as a consistency check.
    """
    val2 = val
    val = posify(val, u2, v2, w2, var2=_vars.var_g1, var3=_vars.var_g2, msg=True, optimize=False)
    if expand(val - val2) != 0:
        raise Exception(f"{val=} {val2=} {u2=} {v2=} {w2=}")

    return val


@cache
def schubmult_generic_partial_posify(u2: Permutation, v2: Permutation):
    """Manifestly positive expansion of ``S_{u2}(x, var_g1) * S_{v2}(x, var_g2)``.

    Returns ``{w2: coeff}`` where each ``coeff`` is the positive representation
    (via ``posify_generic_partial``) of the corresponding
    ``schubmult_double_pair_generic_alt`` coefficient.
    """
    return {w2: posify_generic_partial(val, u2, v2, w2) for w2, val in schubmult_double_pair_generic_alt(u2, v2).items()}


def forwardcoeff(u: Permutation, v: Permutation, perm: Permutation, var2: Alphabet, var3: Alphabet):
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


def dualcoeff(u: Permutation, v: Permutation, perm: Permutation, var2: Alphabet, var3: Alphabet):
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
