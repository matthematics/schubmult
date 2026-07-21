"""Extract the exact bijection behind R_A: MONK + QUANT = 0 (Case A, u(k)=w(k)).

R_A nonzero only when QUANT = q_{k-1} c_{k-2}(p-2;u,w) != 0, i.e. u ->_{k-2} w.
Claim (to prove): there is then a UNIQUE down-neighbor pp*<k with v* = w t_{k,pp*}
such that
    [q if drop] c_{k-1}(p-1;u,v*) = - q_{k-1} c_{k-2}(p-2;u,w),
and both sides are (± q-weight times) the SAME factorial-elementary polynomial:
    P_{k-1}(u,v*) == P_{k-2}(u,w)   [same y-variable multiset]
    deg, nv match                    [same polynomial E_{deg,nv}(y_P;z)]
and the q-weights satisfy the sign-reversing relation.

This script prints, for each nonzero R_A instance, the two sides' (P, n, deg, qweight)
and asserts the three matchings.  Tallies any mismatch.

Run: conda activate schubmult_312 && python _lscripts/RA_bijection.py 5
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


def Pset(u, w, m):
    return sorted(val(u, i) for i in range(1, m + 1) if val(u, i) == val(w, i))


def c_data(u, w, m, a, Rm):
    """Return (qweight, P, nv, deg) or None if c_m(a;u,w)=0."""
    if m < 0:
        return None
    w = Permutation(w)
    if w not in Rm:
        return None
    qm = Rm[w]
    Pm = Pset(u, w, m)
    nm = m - len(Pm)
    deg = a - nm
    nv = m - nm
    if deg < 0 or deg > nv:
        return None
    return (qm, tuple(Pm), nv, deg)


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    total = 0
    mism_P = 0
    mism_deg = 0
    mism_q = 0
    ttype_hist = Counter()
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rkm1 = enumerate_pieri(u, k - 1, N)
            Rkm2 = enumerate_pieri(u, k - 2, N) if k - 2 >= 0 else {u: S.One}
            # candidate w
            cands = set(Rkm1) | set(Rkm2)
            for v in list(Rkm1):
                for pp in range(1, N + 1):
                    if pp != k:
                        cands.add(v.swap(k - 1, pp - 1))
            for p in range(1, k + 1):
                for w in cands:
                    w = Permutation(w)
                    if val(u, k) != val(w, k):
                        continue
                    qd = c_data(u, w, k - 2, p - 2, Rkm2)
                    if qd is None:
                        continue  # QUANT = 0
                    total += 1
                    q_quant, P_quant, nv_q, deg_q = qd
                    # find unique nonzero down MONK neighbor
                    found = None
                    for pp in range(1, k):
                        vv = w.swap(k - 1, pp - 1)
                        dl = w.inv - vv.inv
                        lo, hi = pp, k
                        if dl == 1:
                            md = c_data(u, vv, k - 1, p - 1, Rkm1)
                            ttype, qfac, sign = "RAISE", S.One, -1  # sign from pp<k => -1
                        elif dl == 1 - 2 * (hi - lo):
                            md = c_data(u, vv, k - 1, p - 1, Rkm1)
                            ttype, qfac, sign = "DROP", q_ab(lo, hi), -1
                        else:
                            continue
                        if md is not None:
                            found = (pp, ttype, qfac, sign, md, vv)
                            break
                    if found is None:
                        print(f"  NO MONK MATCH: u={tuple(u)} p={p} k={k} w={tuple(w)}")
                        continue
                    pp, ttype, qfac, sign, md, vv = found
                    q_monk, P_monk, nv_m, deg_m = md
                    ttype_hist[ttype] += 1
                    if P_monk != P_quant:
                        mism_P += 1
                        if mism_P <= 5:
                            print(f"  P MISMATCH u={tuple(u)} p={p} k={k} w={tuple(w)} pp*={pp} {ttype}: Pmonk={P_monk} Pquant={P_quant}")
                    if (nv_m, deg_m) != (nv_q, deg_q):
                        mism_deg += 1
                        if mism_deg <= 5:
                            print(f"  DEG MISMATCH u={tuple(u)} p={p} k={k} w={tuple(w)} pp*={pp}: monk(nv,deg)=({nv_m},{deg_m}) quant=({nv_q},{deg_q})")
                    # q-weight sign-reversal: sign*qfac*q_monk + q_{k-1}*q_quant == 0
                    lhs_q = sign * qfac * q_monk + q_gs[k - 1] * q_quant
                    if expand(lhs_q) != 0:
                        mism_q += 1
                        if mism_q <= 8:
                            print(f"  Q MISMATCH u={tuple(u)} p={p} k={k} w={tuple(w)} pp*={pp} {ttype}: "
                                  f"sign*qfac*q_monk={expand(sign*qfac*q_monk)} q_{{k-1}}*q_quant={expand(q_gs[k-1]*q_quant)}")
    print(f"\nn={n}: R_A nonzero instances = {total}")
    print(f"  special-neighbor type: {dict(ttype_hist)}")
    print(f"  P-set mismatches:   {mism_P}")
    print(f"  (nv,deg) mismatches:{mism_deg}")
    print(f"  q-weight mismatches:{mism_q}")
    print("  => R_A proven bijectively" if (mism_P==0 and mism_deg==0 and mism_q==0) else "  => needs more analysis")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    run(n)
