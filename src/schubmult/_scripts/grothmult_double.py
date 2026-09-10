import sys

from schubmult import GeneratingSet, Permutation, uncode
from schubmult.abc import beta
from schubmult.mult.groth_double import grothmult_double, mult_poly_groth_double
from schubmult.symbolic import S, sstr, sympify, sympify_sympy, expand
from schubmult.utils.argparse import schub_argparse
from schubmult.utils.logging import get_logger

logger = get_logger(__name__)

_posify_cache = {}

# graceful failure ceiling: candidate enumeration and the LP grow combinatorially,
# and unbounded growth has OOM-killed the whole machine before
_MAX_CANDIDATES = 1_000_000


def _limit_memory():
    """Cap this process's address space so runaway cases die with MemoryError
    instead of taking down the host (WSL is especially fragile under OOM)."""
    import resource

    try:
        import psutil

        cap = max(int(psutil.virtual_memory().available * 0.7), 1024**3)
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        if soft == resource.RLIM_INFINITY or soft > cap:
            resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    except Exception:  # noqa: BLE001
        logger.debug("could not set memory limit", exc_info=True)


def _solver(msg):
    import pulp as pu

    if pu.HiGHS_CMD().available():
        return pu.HiGHS_CMD(msg=msg)
    return pu.PULP_CBC_CMD(msg=msg)


def _denominator_budget(den, var2):
    """Atom multiplicities ``{y_symbol: exp}`` of a denominator ``prod (1 + beta*y_i)**e``.

    Keys are symengine symbols (consistent with GeneratingSet identity).
    """
    from sympy import factor_list

    bs = sympify_sympy(beta)
    budget = {}
    _, factors = factor_list(sympify_sympy(den))
    for fac, mult in factors:
        syms = [s for s in fac.free_symbols if s != bs]
        if len(syms) != 1 or var2.index(syms[0]) == -1:
            raise ValueError(f"unexpected denominator factor {fac}")
        key = sympify(syms[0])
        budget[key] = budget.get(key, 0) + int(mult)
    return budget


def groth_posify(val, var2, var3, msg):
    """Positive FGL form of a coefficient, at ``beta = 1``.

    Basis: all products of ``z_i - y_j``, ``(1 + y_j)^{-1}``, ``(1 + z_i)^{-1}``.
    Substituting ``u_j = 1 + y_j``, ``v_i = 1 + z_i`` (so ``y = u - 1``,
    ``z = v - 1``) turns the value into a Laurent polynomial in ``u, v``; clearing
    the monomial denominator makes everything polynomial.  Candidates are
    difference products times leftover clearing monomials, expanded once, with
    ``as_coefficients_dict`` terms used as opaque basis vectors for an integer
    LP.  ``beta`` is restored as ``beta**(#diffs - d)`` with the Laurent atoms
    ``(1 + beta*y)`` degree 0 by construction.
    """
    import itertools

    import pulp as pu
    from symengine import Symbol as SySymbol

    from schubmult.symbolic import efficient_subs

    try:
        return int(expand(val))
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

    uu = {ys: SySymbol(f"u_{var2.index(ys)}") for ys in ysyms}
    vv = {zs: SySymbol(f"v_{var3.index(zs)}") for zs in zsyms}
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
            key = 0
            deg = 0
            for g, p in mono.as_powers_dict().items():
                if g in gshift:
                    key += int(p) << gshift[g]
                    deg += int(p)
                elif g != S.One:
                    raise ValueError(f"unexpected generator {g}")
            vec[key] = vec.get(key, 0) + int(c)
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
    # candidates are homogeneous of total degree #pairs + clearing-shift degree,
    # so degrees absent from the target can be pruned losslessly
    target_degs = frozenset(deg0.values())

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

    def get_level(m):
        if m in state["levels"]:
            return state["levels"][m]
        prods = prods_state["prods"]
        pair_vecs = prods_state["pair_vecs"]
        shifts = [(shift, usage) for shift, usage, sdeg in state["usage_shifts"] if m + sdeg in target_degs]
        lev = []
        if shifts:
            while len(prods) <= m:
                cur = []
                for vec, combo, start in prods[-1]:
                    for idx in range(start, len(pair_vecs)):
                        pv, pair = pair_vecs[idx]
                        cur.append((dict_mul(vec, pv), (*combo, pair), idx))
                prods.append(cur)
            for vec, combo, _ in prods[m]:
                for shift, usage in shifts:
                    lev.append(({k + shift: c for k, c in vec.items()}, combo, usage, m))
        state["levels"][m] = lev
        return lev

    def solve(candidates):
        vrs = [pu.LpVariable(name=f"a{i}", lowBound=0, cat="Integer") for i in range(len(candidates))]
        lp_prob = pu.LpProblem("Problem", pu.LpMinimize)
        lp_prob += 0
        eqs = {}
        for i, (svec, _, _, _) in enumerate(candidates):
            for k, c in svec.items():
                eqs.setdefault(k, {})[vrs[i]] = c
        for k in set(eqs) | set(vec0):
            lp_prob += pu.LpAffineExpression(eqs.get(k, {})) == vec0.get(k, 0)
        try:
            status = lp_prob.solve(_solver(msg))
        except KeyboardInterrupt:
            import psutil

            current_process = psutil.Process()
            for child in current_process.children(recursive=True):
                child_process = psutil.Process(child.pid)
                child_process.terminate()
                child_process.kill()
            raise
        return status, vrs

    # escalate the factor-count ceiling: small LPs solve fast and usually suffice
    candidates = []
    status = None
    for m in range(max(d, 0), max_m + 1):
        lev = get_level(m)
        if not lev:
            continue
        candidates = [*candidates, *lev]
        if len(candidates) > _MAX_CANDIDATES:
            raise ValueError(f"candidate set too large ({len(candidates)}) for {val}")
        print(f"  solving level m={m}: {len(candidates)} candidates", file=sys.stderr)
        status, vrs = solve(candidates)
        if pu.LpStatus[status] == "Optimal":
            break
    if status is None or pu.LpStatus[status] != "Optimal":
        raise ValueError(f"no positive representation found for {val}")

    result = S.Zero
    for i, (_, combo, usage, m) in enumerate(candidates):
        x = vrs[i].value()
        # round, don't truncate: solvers return near-integers like 0.9999999999996
        xi = 0 if x is None else round(x)
        if xi != 0:
            term = S.One
            for ys, zs in combo:
                term = term * (sympify(zs) - sympify(ys))
            for g, a in usage.items():
                if a:
                    term = term / back[g] ** a
            result += xi * bs ** (m - d) * term
    # exact symbolic verification: structural numerator of the difference
    diff_num, _ = sympify_sympy(result - e).as_numer_denom()
    if expand(sympify(diff_num)) != S.Zero:
        raise ValueError(f"positive representation check failed for {val}: got {result}")
    return result


def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        var2 = GeneratingSet("y")
        var3 = GeneratingSet("z")
        sys.setrecursionlimit(1000000)

        args, formatter = schub_argparse(
            "grothmult_double",
            "Compute coefficients of products of double Grothendieck polynomials in the same or different sets of coefficient variables",
            argv=argv[1:],
            yz=True,
            coprod=False,
        )

        if args.display_positive and args.same:
            print("--display-positive is only supported with --mixed-var for grothmult_double")
            return 1

        if args.display_positive:
            _limit_memory()

        mult = args.mult
        mulstring = args.mulstring

        perms = args.perms

        ascode = args.ascode
        Permutation.print_as_code = ascode
        same = args.same
        pr = args.pr

        if same:
            var3 = var2

        if ascode:
            perms = [uncode(perm) for perm in perms]
        else:
            for i in range(len(perms)):
                if len(perms[i]) < 2 and (len(perms[i]) == 0 or perms[i][0] == 1):
                    perms[i] = Permutation([])
                perms[i] = Permutation(perms[i])

        coeff_dict = {perms[0]: 1}

        for perm in perms[1:]:
            coeff_dict = grothmult_double(coeff_dict, perm, var2, var3)

        if mult:
            mul_exp = sympify(mulstring)
            coeff_dict = mult_poly_groth_double(coeff_dict, mul_exp, var2, var3)

        # sort/filter up front so posified coefficients can stream out as each
        # one finishes (the LPs can take a long time)
        coeff_perms = [perm for perm, val in coeff_dict.items() if expand(val) != 0]
        coeff_perms.sort(key=lambda x: (-abs(perms[0].inv + perms[1].inv - x.inv), *x))
        width = max([len(sstr(perm)) for perm in coeff_perms]) if coeff_perms else 0

        raw_result_dict = {}
        for i, perm in enumerate(coeff_perms):
            val = coeff_dict[perm]
            if args.display_positive:
                print(f"posify {i + 1}/{len(coeff_perms)}: {sstr(perm)}", file=sys.stderr)
                try:
                    # groth_posify verifies its own output exactly
                    val = groth_posify(val, var2, var3, args.msg)
                except Exception:
                    import traceback

                    traceback.print_exc()
                    print(f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=}")
                    return 1
            raw_result_dict[perm] = val
            if pr and formatter:
                print(f"{sstr(perm)!s:>{width}}  {formatter(val)}", flush=True)

        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))
