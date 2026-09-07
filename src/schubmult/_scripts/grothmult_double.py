import sys

from schubmult import GeneratingSet, Permutation, uncode
from schubmult.abc import beta
from schubmult.mult.groth_double import grothmult_double, mult_poly_groth_double
from schubmult.symbolic import S, sstr, sympify, sympify_sympy, expand
from schubmult.utils.argparse import schub_argparse
from schubmult.utils.logging import get_logger

logger = get_logger(__name__)


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

    vec0 = {k: int(c) for k, c in num_l.as_coefficients_dict().items() if c != 0}
    max_m = max((sum(int(p) for p in mono.as_powers_dict().values()) for mono in vec0), default=0)

    # candidates: product of m difference atoms times a leftover clearing
    # monomial (unused inverse depth); expand and keep the coefficient dict
    pairs = [(ys, zs) for ys in ysyms for zs in zsyms]
    cap_syms = list(caps)
    candidates = []
    for m in range(max(d, 0), max_m + 1):
        for combo in itertools.combinations_with_replacement(pairs, m):
            prd = S.One
            for ys, zs in combo:
                prd = prd * (vv[zs] - uu[ys])
            for usage in itertools.product(*[range(caps[g] + 1) for g in cap_syms]):
                cand = prd
                for g, a in zip(cap_syms, usage):
                    cand = cand * g ** (caps[g] - a)
                svec = {k: int(c) for k, c in expand(cand).as_coefficients_dict().items() if c != 0}
                candidates.append((svec, combo, dict(zip(cap_syms, usage)), m))

    vrs = [pu.LpVariable(name=f"a{i}", lowBound=0, cat="Integer") for i in range(len(candidates))]
    lp_prob = pu.LpProblem("Problem", pu.LpMinimize)
    lp_prob += 0
    eqs = {}
    for i, (svec, _, _, _) in enumerate(candidates):
        for k, c in svec.items():
            eqs.setdefault(k, []).append(c * vrs[i])
    for k in set(eqs) | set(vec0):
        lp_prob += pu.lpSum(eqs.get(k, [])) == vec0.get(k, 0)
    try:
        solver = pu.PULP_CBC_CMD(msg=msg)
        status = lp_prob.solve(solver)
    except KeyboardInterrupt:
        import psutil

        current_process = psutil.Process()
        for child in current_process.children(recursive=True):
            child_process = psutil.Process(child.pid)
            child_process.terminate()
            child_process.kill()
        raise
    if pu.LpStatus[status] != "Optimal":
        raise ValueError(f"no positive representation found for {val}")

    result = S.Zero
    for i, (_, combo, usage, m) in enumerate(candidates):
        x = vrs[i].value()
        if x is not None and int(x) != 0:
            term = S.One
            for ys, zs in combo:
                term = term * (sympify(zs) - sympify(ys))
            for g, a in usage.items():
                if a:
                    term = term / back[g] ** a
            result += int(x) * bs ** (m - d) * term
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

        if args.display_positive:
            new_dict = {}
            for perm, val in coeff_dict.items():
                if expand(val) == 0:
                    continue
                try:
                    # groth_posify verifies its own output exactly
                    pos_val = groth_posify(val, var2, var3, args.msg)
                except Exception:
                    import traceback

                    traceback.print_exc()
                    print(f"error; write to schubmult@gmail.com with the case {perms=} {perm=} {val=}")
                    return 1
                new_dict[perm] = pos_val
            coeff_dict = new_dict

        raw_result_dict = {}
        if pr or formatter is None:
            width = max([len(sstr(perm)) for perm in coeff_dict]) if coeff_dict else 0
            coeff_perms = list(coeff_dict.keys())
            coeff_perms.sort(key=lambda x: (-abs(perms[0].inv + perms[1].inv - x.inv), *x))

            for perm in coeff_perms:
                val = coeff_dict[perm]
                if expand(val) != 0:
                    raw_result_dict[perm] = val
                    if formatter:
                        print(f"{sstr(perm)!s:>{width}}  {formatter(val)}")

        if formatter is None:
            return raw_result_dict
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))
