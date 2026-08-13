"""rho<=1 via log-concavity / LGV, NOT via subposet restriction.

Setup: pad all of alpha, beta, gamma=alpha+beta to a common ambient S_M by fixed points.
Lambda(sigma_c) = M! / prod_{i=1}^M h_i(c),  h_i(c) = down-set size of node i in T_c
              = 1 + #{ j<i : sigma_c(j) < sigma_c(i) }  = sigma_c(i) - c_i  (with padding).
Padding by fixed points multiplies numerator by the new factors and denominator by the
same (h=full for a fixed point at the end), leaving Lambda unchanged -> ambient-independent.

Then
   rho = Lambda(gamma)/(Lambda(alpha)Lambda(beta))
       = [M!/prod g_i] / ([M!/prod a_i][M!/prod b_i])
       = (1/M!) * prod_i ( a_i b_i / g_i ).
So rho<=1  <=>  prod_i a_i b_i  <=  M! * prod_i g_i.           (HOOK INEQUALITY)

TESTS:
 (T1) verify the identity rho = (1/M!) prod a_i b_i/g_i exactly.
 (T2) test the hook inequality prod a_i b_i <= M! prod g_i.
 (T3) LOG-CONCAVITY probe: is  Lambda(c)^2 >= Lambda(c-e_i)Lambda(c+e_i)  (discrete concavity
      of log Lambda along each coordinate)?  This is the standard route to submultiplicativity.
 (T4) does submultiplicativity follow from (T3)+monotone marginals along a MONOTONE path
      c: 0 -> alpha -> alpha+beta adding beta's boxes in NONINCREASING row order?
"""
from math import factorial
from itertools import product as iproduct
from fractions import Fraction
from schubmult import uncode


def sigma_oneline(code, M):
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    if not code:
        return list(range(1, M+1))
    ol = list(uncode(code))
    return ol + list(range(len(ol)+1, M+1))


def hooks(code, M):
    ol = sigma_oneline(code, M)
    return [ol[i] - (code[i] if i < len(code) else 0) for i in range(M)]


def Lam(code, M):
    h = hooks(code, M)
    p = 1
    for x in h:
        p *= x
    return factorial(M) // p


def is_part(t):
    return all(t[i] >= t[i+1] for i in range(len(t)-1)) and all(x >= 0 for x in t)


def main():
    M = 7
    maxpart, maxlen = 4, 4
    parts = [t for t in iproduct(range(maxpart+1), repeat=maxlen) if is_part(t)]

    # T1 + T2
    id_fail = hook_fail = 0
    worst_rho = None
    n_pairs = 0
    for a in parts:
        for b in parts:
            L = maxlen
            g = tuple(a[i]+b[i] for i in range(L))
            if max(g) > M-1:
                continue
            ha = hooks(a, M); hb = hooks(b, M); hg = hooks(g, M)
            pa = 1; pb = 1; pg = 1
            for i in range(M):
                pa *= ha[i]; pb *= hb[i]; pg *= hg[i]
            rho = Fraction(Lam(g, M), Lam(a, M)*Lam(b, M))
            rho_hooks = Fraction(pa*pb, factorial(M)*pg)
            if rho != rho_hooks:
                id_fail += 1
            if pa*pb > factorial(M)*pg:
                hook_fail += 1
            if worst_rho is None or rho > worst_rho[0]:
                worst_rho = (rho, a, b)
            n_pairs += 1
    print(f"T1 identity rho=(1/M!)prod a_i b_i/g_i : {id_fail} failures / {n_pairs}")
    print(f"T2 hook ineq prod a_i b_i <= M! prod g_i: {hook_fail} failures / {n_pairs}; sup rho={float(worst_rho[0]):.4f} at a={worst_rho[1]} b={worst_rho[2]}")

    # T3 discrete log-concavity along each coordinate
    lc_fail = 0; lc_checks = 0; lc_worst = None
    for c in parts:
        for i in range(maxlen):
            cm = list(c); cm[i] -= 1
            cp = list(c); cp[i] += 1
            if not is_part(cm) or not is_part(cp):
                continue
            if max(cp) > M-1:
                continue
            lhs = Lam(c, M)**2
            rhs = Lam(tuple(cm), M) * Lam(tuple(cp), M)
            lc_checks += 1
            if lhs < rhs:
                lc_fail += 1
                r = Fraction(rhs, lhs)
                if lc_worst is None or r > lc_worst[0]:
                    lc_worst = (r, c, i)
    print(f"T3 coordinate log-concavity Lam(c)^2>=Lam(c-e_i)Lam(c+e_i): {lc_fail} failures / {lc_checks}")
    if lc_worst:
        print(f"     worst: c={lc_worst[1]} i={lc_worst[2]} ratio={float(lc_worst[0]):.4f}")

    # T4 monotone-path submultiplicativity: add beta boxes top row first (noninc) vs any order
    # Test: along path 0 -> alpha -> alpha+beta adding beta boxes, does each marginal
    #   R_i(current) <= R_i(current - alpha-part) ... just test the end-to-end product telescopes <=1
    # (already implied by T2; T4 instead checks a SUFFICIENT local condition:
    #  adding a box never increases Lambda by more than the standalone factor)
    # marginal on empty base for one box in row i:
    def Ri(code, i, M):
        c2 = list(code); c2[i] += 1
        if not is_part(c2) or max(c2) > M-1:
            return None
        return Fraction(Lam(tuple(c2), M), Lam(code, M))
    # SUFFICIENT: R_i(alpha + partial_beta) <= R_i(partial_beta) when boxes added in
    # noninc row order AND partial_beta has been built the same way. Test that specific order.
    t4_fail = 0; t4_checks = 0; t4_worst=None
    for a in parts:
        for b in parts:
            L = maxlen
            if max(a[i]+b[i] for i in range(L)) > M-1:
                continue
            # build beta by rows in NONINCREASING row index order = row 0 fully, then row1...
            base_b = [0]*L
            base_ab = list(a)
            ok = True
            for i in range(L):
                for _ in range(b[i]):
                    ri_heavy = Ri(tuple(base_ab), i, M)
                    ri_light = Ri(tuple(base_b), i, M)
                    if ri_heavy is None or ri_light is None:
                        ok = False; break
                    t4_checks += 1
                    if ri_heavy > ri_light:
                        t4_fail += 1
                        d = ri_heavy - ri_light
                        if t4_worst is None or d > t4_worst[0]:
                            t4_worst = (d, a, b, i, ri_heavy, ri_light)
                    base_ab[i] += 1
                    base_b[i] += 1
                if not ok:
                    break
    print(f"T4 row-order marginal antitone (heavy<=light along fill): {t4_fail} failures / {t4_checks}")
    if t4_worst:
        d,a,b,i,rh,rl = t4_worst
        print(f"     worst: a={a} b={b} row={i} R_heavy={float(rh):.3f} > R_light={float(rl):.3f}")


if __name__ == "__main__":
    main()
