"""Nail the DETERMINISTIC rule for pp* in the R_A bijection, to make it constructive/provable.

Hypotheses to test (Case A, u(k)=w(k), QUANT!=0, unique down-neighbor pp*, v=w t_{pp*,k}):
  H1: pp* = k-1 - (# of consecutive i just below k with u(i)=w(i)), i.e. skip fixed points.
      Equivalently pp* = largest j<k with u(j)!=w(j), OR ... test a few variants.
  H2: v is obtained from w by a cyclic shift of the block [pp*..k].
  H3: relation between v and w on the block.

Prints, per instance: k, w restricted to [1..k], fixed-point pattern on [1..k-1], pp*, v[1..k].
Then tests candidate closed formulas and reports match rate.

Run: conda activate schubmult_312 && python _lscripts/ppstar_rule.py 5
"""

import sys
from collections import Counter

from schubmult import *  # noqa: F401,F403
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, expand_func, prod
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import FactorialElemSym

y = GeneratingSet("y")
z = GeneratingSet("z")
q_gs = GeneratingSet("q")


def q_ab(a, b):
    return prod([q_gs[s] for s in range(a, b)])


def val(perm, pos):
    return perm[pos - 1]


def enumerate_pieri(u, k, N):
    u = Permutation(u)
    results = {u: S.One}
    stack = [(u, frozenset(), N + 1, S.One, u.inv)]
    seen = set()
    while stack:
        perm, used_a, last_b, qw, clen = stack.pop()
        key = (perm, used_a, last_b)
        if key in seen:
            continue
        seen.add(key)
        for a in range(1, k + 1):
            if a in used_a:
                continue
            for b in range(k + 1, N + 1):
                if b > last_b:
                    continue
                nperm = perm.swap(a - 1, b - 1)
                d = nperm.inv - clen
                if d == 1:
                    nqw = qw
                elif d == -2 * (b - a) + 1:
                    nqw = qw * q_ab(a, b)
                else:
                    continue
                results.setdefault(nperm, nqw)
                stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results


def P_set(u, w, m):
    return sorted(val(u, i) for i in range(1, m + 1) if val(u, i) == val(w, i))


def c_term(u, w, m, a, Rm):
    if m < 0:
        return S.Zero
    w = Permutation(w)
    if w not in Rm:
        return S.Zero
    qm = Rm[w]
    Pm = P_set(u, w, m)
    nm = m - len(Pm)
    deg = a - nm
    nv = m - nm
    if deg < 0 or deg > nv:
        return S.Zero
    if nv == 0:
        return qm if deg == 0 else S.Zero
    yvars = [y[v] for v in Pm]
    ncoeff = nv + 1 - deg
    zvars = [z[i] for i in range(1, ncoeff + 1)]
    return qm * expand_func(FactorialElemSym(deg, nv, yvars, zvars))


def candidate_ppstar(u, w, k):
    # H1: largest j<k with u(j)!=w(j); if none, ... fall back
    js = [j for j in range(1, k) if val(u, j) != val(w, j)]
    h1 = max(js) if js else None
    return h1


def run(n, show=25):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    h1_ok = 0
    h1_bad = 0
    total = 0
    shown = 0
    examples = []
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rkm1 = enumerate_pieri(u, k - 1, N)
            Rkm2 = enumerate_pieri(u, k - 2, N) if k - 2 >= 0 else {u: S.One}
            cands = set(Rkm1) | set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1, N + 1):
                    if pp != k:
                        cands.add(v0.swap(k - 1, pp - 1))
            for p in range(1, k + 1):
                for w in cands:
                    w = Permutation(w)
                    if val(u, k) != val(w, k):
                        continue
                    QUANT = q_gs[k - 1] * c_term(u, w, k - 2, p - 2, Rkm2)
                    if expand(QUANT) == 0:
                        continue
                    downs = []
                    for pp in range(1, k):
                        vv = w.swap(k - 1, pp - 1)
                        dl = w.inv - vv.inv
                        lo, hi = pp, k
                        if dl == 1 or dl == 1 - 2 * (hi - lo):
                            cval = c_term(u, vv, k - 1, p - 1, Rkm1)
                            if expand(cval) != 0:
                                downs.append(pp)
                    if len(downs) != 1:
                        continue
                    ppstar = downs[0]
                    total += 1
                    h1 = candidate_ppstar(u, w, k)
                    if h1 == ppstar:
                        h1_ok += 1
                    else:
                        h1_bad += 1
                        if len(examples) < show:
                            fixed = [1 if val(u, i) == val(w, i) else 0 for i in range(1, k)]
                            examples.append((tuple(u), p, k, tuple(w)[:k + 1], ppstar, h1, fixed))
    print(f"n={n}: R_A unique-neighbor instances={total}")
    print(f"  H1 (largest j<k with u(j)!=w(j)) == pp*:  {h1_ok}/{total}   bad={h1_bad}")
    for ex in examples:
        print(f"    u={ex[0]} p={ex[1]} k={ex[2]} w[:k+1]={ex[3]} pp*={ex[4]} H1={ex[5]} fixed<k={ex[6]}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    run(n)
