#!/usr/bin/env python3
"""Unconditional bounds for schubmult.

Quantities (v in S_n, b = mdom v, theta = theta(v^{-1})):
    D(v) = prod_i (theta_i + 1)        L = #{i : theta_i != 0}      D <= n^L
    Lambda(w) = |[e,w]_L|
    N(u,v)   = # nonzero terms of Sch_u(x;y) Sch_v(x;z)
             = | union_{v'' <=_L v} supp( Sch_u Sch_{v''} ) |   (double support thm)

Tests:
  T1  N(u,v) from the double product == the union-of-single-supports formula
  T2  N(e,v) == Lambda(v)                   [=> u=e is the extremal instance]
  T3  Lambda(mdom v) <= D(v)                [the paper's commented domhullbound]
  T4  N(u, mdom v) <= D(v) * N(u,v)         [sharpening of N^up <= n^L N]
  T5  the separating family a_k = s_2 s_4 ... s_{2k} in S_{2k+1}
"""

import sys
from collections import defaultdict

sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")

from schubmult import DSx, Sx  # noqa: E402
from schubmult.combinatorics.permutation import Permutation, uncode  # noqa: E402
from schubmult.symbolic import S, expand, sympify  # noqa: E402


def D_of(v):
    d = 1
    L = 0
    for t in (~v).theta():
        if t:
            d *= t + 1
            L += 1
    return d, L


def covers_down(w):
    wi = ~w
    return [~(wi.swap(d, d + 1)) for d in wi.zero_indexed_descents()]


def below(w):
    seen = {w}
    fr = [w]
    while fr:
        nxt = []
        for x in fr:
            for y in covers_down(x):
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
        fr = nxt
    return seen


_ssupp = {}


def single_supp(u, v):
    if (u, v) not in _ssupp:
        _ssupp[(u, v)] = frozenset((Sx(u) * Sx(v)).keys())
    return _ssupp[(u, v)]


def N_formula(u, v, bel=None):
    bel = bel if bel is not None else below(v)
    out = set()
    for vpp in bel:
        out |= single_supp(u, vpp)
    return len(out)


def N_direct(u, v):
    prod = DSx(u) * DSx(v, "z")
    return sum(1 for c in prod.values() if expand(sympify(c)) != S.Zero)


def t1_t2(n):
    perms = list(Permutation.all_permutations(n))
    bad1 = bad2 = 0
    for v in perms:
        bel = below(v)
        if N_formula(Permutation([]), v, bel) != len(bel):
            bad2 += 1
        for u in perms:
            if N_formula(u, v, bel) != N_direct(u, v):
                bad1 += 1
                if bad1 <= 2:
                    print(f"      T1 FAIL u={u} v={v}")
    print(f"   T1 (n={n}) formula vs double product: {bad1} mismatches")
    print(f"   T2 (n={n}) N(e,v) == Lambda(v):       {bad2} mismatches")


def t3(hi):
    print("\n   T3  Lambda(mdom v) <= D(v) = prod(theta_i+1)")
    for n in range(2, hi + 1):
        worst = 0.0
        eqs = 0
        tot = 0
        seen_b = {}
        for v in Permutation.all_permutations(n):
            b = v.mul_dominant()
            if b not in seen_b:
                seen_b[b] = len(below(b))
            lb = seen_b[b]
            d, _L = D_of(v)
            tot += 1
            if lb == d:
                eqs += 1
            worst = max(worst, lb / d)
        flag = "OK" if worst <= 1 else "*** VIOLATED ***"
        print(f"      n={n}: max Lambda(b)/D = {worst:.4f}  {flag}   (equality in {eqs}/{tot})")


def t4(hi):
    print("\n   T4  N(u, mdom v) <= D(v) * N(u,v)   [vs the paper's n^L]")
    for n in range(2, hi + 1):
        perms = list(Permutation.all_permutations(n))
        bel = {v: below(v) for v in perms}
        worstD = worstL = 0.0
        for v in perms:
            b = v.mul_dominant()
            d, L = D_of(v)
            for u in perms:
                r = N_formula(u, b, bel[b]) / N_formula(u, v, bel[v])
                worstD = max(worstD, r / d)
                worstL = max(worstL, r / n**L)
        fD = "OK" if worstD <= 1 else "*** VIOLATED ***"
        fL = "OK" if worstL <= 1 else "*** VIOLATED ***"
        print(f"      n={n}: max ratio/D = {worstD:.4f} {fD}      max ratio/n^L = {worstL:.4f} {fL}")


def t5(kmax):
    print("\n   T5  separating family  a_k = s_2 s_4 ... s_{2k}  in  S_{2k+1}")
    print("       k   n    L    D(v)    N(e,v)=Lambda(v)   Lambda(b)   log N(work)/log N(out)")
    import math

    for k in range(1, kmax + 1):
        n = 2 * k + 1
        code = [0, 1] * k
        v = uncode(code)
        d, L = D_of(v)
        lv = len(below(v))
        lb = len(below(v.mul_dominant()))
        ratio = math.log(lb) / math.log(lv) if lv > 1 else float("inf")
        print(f"       {k}  {n:>3}  {L:>3}  {d:>7}   {lv:>14}   {lb:>9}   {ratio:>8.4f}")


if __name__ == "__main__":
    t1_t2(4)
    t3(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
    t4(5)
    t5(5)
