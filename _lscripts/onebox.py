"""Test the one-box marginal lemma underlying rho<=1.

Lambda(sigma_lambda) via closed form Lambda(u)=n!/prod_i(u(i)-code_i(u)),
computed on the dominant permutation with code = partition lambda.
Ambient-independent (padding by fixed points leaves quotient unchanged).

Define R_i(gamma) = Lambda(sigma_{gamma+e_i}) / Lambda(sigma_gamma).

CLAIM (antitone marginal / diminishing returns):
   gamma' <= gamma coordinatewise  =>  R_i(gamma) <= R_i(gamma')
for all rows i such that gamma+e_i and gamma'+e_i are both partitions.

Also directly re-verify:  rho = Lambda(a+b)/(Lambda(a)Lambda(b)) <= 1.
"""
import sys
from math import factorial
from itertools import product as iproduct
from fractions import Fraction
from schubmult import Permutation, uncode


def Lam(code):
    """Lambda(sigma) for dominant sigma with given code (a partition, trailing zeros ok)."""
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    L = len(code)
    if L == 0:
        return 1
    n = L + (code[0] if code else 0) + 1
    perm = uncode(code)
    ol = list(perm)  # one-line
    ol = ol + list(range(len(ol) + 1, n + 1))
    c = list(perm.code) + [0] * (n - len(perm.code))
    denom = 1
    for i in range(n):
        denom *= (ol[i] - c[i])
    return factorial(n) // denom


def is_partition(t):
    return all(t[i] >= t[i + 1] for i in range(len(t) - 1)) and all(x >= 0 for x in t)


def partitions_upto(maxpart, maxlen):
    out = []
    for tup in iproduct(range(maxpart + 1), repeat=maxlen):
        if is_partition(tup):
            out.append(tup)
    return out


def main():
    maxpart, maxlen = 4, 4
    parts = partitions_upto(maxpart, maxlen)
    lam = {p: Lam(p) for p in parts}

    # sanity: rectangles give binomials
    for (c, k) in [(2, 2), (3, 3), (2, 3), (1, 4)]:
        rect = tuple([c] * k)
        from math import comb
        assert Lam(rect) == comb(c + k, k), (rect, Lam(rect), comb(c + k, k))
    print("rectangle closed-form check passed")

    # rho <= 1
    bad_rho = 0
    worst_rho = None
    P = parts
    cnt = 0
    for a in P:
        for b in P:
            L = maxlen
            g = tuple(a[i] + b[i] for i in range(L))
            if max(g) > maxpart:  # keep g in-table by recomputing directly
                lg = Lam(g)
            else:
                lg = lam[g]
            r = Fraction(lg, lam[a] * lam[b])
            cnt += 1
            if r > 1:
                bad_rho += 1
            if worst_rho is None or r > worst_rho[0]:
                worst_rho = (r, a, b, lg, lam[a], lam[b])
    print(f"rho<=1: {bad_rho} failures over {cnt} pairs; sup rho = {float(worst_rho[0]):.4f} at a={worst_rho[1]} b={worst_rho[2]}")

    # one-box marginal antitone test
    def Ri(gamma, i):
        g2 = list(gamma)
        g2[i] += 1
        if not is_partition(g2):
            return None
        return Fraction(Lam(g2), Lam(gamma))

    bad_anti = 0
    worst_anti = None
    checks = 0
    for gp in parts:
        for g in parts:
            if not all(gp[i] <= g[i] for i in range(maxlen)):
                continue
            for i in range(maxlen):
                rg = Ri(g, i)
                rgp = Ri(gp, i)
                if rg is None or rgp is None:
                    continue
                checks += 1
                # claim R_i(g) <= R_i(g')  (bigger base => smaller marginal)
                if rg > rgp:
                    bad_anti += 1
                    if worst_anti is None or (rg - rgp) > worst_anti[0]:
                        worst_anti = (rg - rgp, g, gp, i, rg, rgp)
    print(f"antitone R_i: {bad_anti} failures over {checks} comparisons")
    if worst_anti:
        d, g, gp, i, rg, rgp = worst_anti
        print(f"  worst violation: i={i} g={g} g'={gp}  R_i(g)={float(rg):.4f} > R_i(g')={float(rgp):.4f}")


if __name__ == "__main__":
    main()
