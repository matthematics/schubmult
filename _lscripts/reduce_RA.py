"""Reduce R_A (Case A: MONK+QUANT=0) to explicit matching identities, and verify them.

In Case A (u(k)=w(k)) with QUANT = q_{k-1} c_{k-2}(p-2;u,w) != 0, we claim there is a UNIQUE
down-neighbor pp* < k with v := w t_{k,pp*} such that c_{k-1}(p-1;u,v) != 0, and that
   qcoeff(pp*) * c_{k-1}(p-1;u,v) == q_{k-1} * c_{k-2}(p-2;u,w)          (*)
(so MONK = -qcoeff*c_{k-1}(p-1;u,v) = -QUANT).  Since c_m(a;u,.)=q_m E_{deg,nv}(y_P;z), (*) holds
IF the following four component identities hold:
   (M1) P_{k-1}(u,v) == P_{k-2}(u,w)            [same y-variable multiset  => same E]
   (M2) n_{k-1}(u,v) == n_{k-2}(u,w) + 1        [=> deg and nv both match]
   (M3) qcoeff(pp*) * q_{k-1}(u,v) == q_{k-1}_generator * q_{k-2}(u,w)   [q-weights match]
And conversely, when QUANT == 0, no down-neighbor contributes (MONK==0).

This script verifies: existence+uniqueness of pp*, (M1),(M2),(M3), and the QUANT==0 => MONK==0
direction, for all Case-A (u,p,k,w).  If all pass, R_A is reduced to a purely combinatorial
bijection statement whose value-identity is automatic.

Run: conda activate schubmult_312 && python _lscripts/reduce_RA.py 5
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


def n_val(u, w, m):
    return m - len(P_set(u, w, m))


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


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    stats = Counter()
    fails = []
    ppstar_offset = Counter()
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rk = enumerate_pieri(u, k, N)
            Rkm1 = enumerate_pieri(u, k - 1, N)
            Rkm2 = enumerate_pieri(u, k - 2, N) if k - 2 >= 0 else {u: S.One}
            cands = set(Rk) | set(Rkm1) | set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1, N + 1):
                    if pp != k:
                        cands.add(v0.swap(k - 1, pp - 1))
            for p in range(1, k + 1):
                for w in cands:
                    w = Permutation(w)
                    if val(u, k) != val(w, k):
                        continue
                    # collect nonzero down MONK neighbors
                    downs = []
                    for pp in range(1, k):
                        vv = w.swap(k - 1, pp - 1)
                        dl = w.inv - vv.inv
                        lo, hi = pp, k
                        if dl == 1:
                            qc = S.One
                        elif dl == 1 - 2 * (hi - lo):
                            qc = q_ab(lo, hi)
                        else:
                            continue
                        cval = c_term(u, vv, k - 1, p - 1, Rkm1)
                        if expand(cval) != 0:
                            downs.append((pp, vv, qc, cval))
                    QUANT = q_gs[k - 1] * c_term(u, w, k - 2, p - 2, Rkm2)
                    quant_zero = expand(QUANT) == 0

                    stats["total"] += 1
                    if quant_zero:
                        # expect MONK (down) == 0, i.e. no down neighbor
                        if len(downs) != 0:
                            stats["FAIL_quantzero_but_down"] += 1
                            fails.append(("quantzero_down", tuple(u), p, k, tuple(w)))
                        continue
                    # QUANT != 0: expect exactly one down neighbor with matching identities
                    stats["quant_nonzero"] += 1
                    if len(downs) != 1:
                        stats["FAIL_not_unique"] += 1
                        fails.append(("not_unique", tuple(u), p, k, tuple(w), len(downs)))
                        continue
                    pp, vv, qc, cval = downs[0]
                    ppstar_offset[pp - (k - 1)] += 1
                    # (M1)
                    P_v = P_set(u, vv, k - 1)
                    P_w = P_set(u, w, k - 2)
                    if P_v != P_w:
                        stats["FAIL_M1"] += 1
                        fails.append(("M1", tuple(u), p, k, tuple(w), P_v, P_w))
                    else:
                        stats["OK_M1"] += 1
                    # (M2)
                    if n_val(u, vv, k - 1) != n_val(u, w, k - 2) + 1:
                        stats["FAIL_M2"] += 1
                        fails.append(("M2", tuple(u), p, k, tuple(w)))
                    else:
                        stats["OK_M2"] += 1
                    # (M3)
                    q_v = Rkm1[vv]
                    q_w = Rkm2[Permutation(w)]
                    if expand(qc * q_v - q_gs[k - 1] * q_w) != 0:
                        stats["FAIL_M3"] += 1
                        fails.append(("M3", tuple(u), p, k, tuple(w), str(qc), str(q_v), str(q_w)))
                    else:
                        stats["OK_M3"] += 1
                    # full value identity (*)
                    if expand(qc * cval - QUANT) != 0:
                        stats["FAIL_value"] += 1
                        fails.append(("value", tuple(u), p, k, tuple(w)))
                    else:
                        stats["OK_value"] += 1
    print(f"n={n}: {dict(stats)}")
    print(f"  pp*-(k-1) offsets: {dict(sorted(ppstar_offset.items()))}")
    for f in fails[:15]:
        print("  FAIL:", f)
    return not any(kk.startswith("FAIL") for kk in stats)


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    ok = run(n)
    print("RESULT:", "R_A reduced & matching identities hold" if ok else "some matching identity FAILS")
    sys.exit(0 if ok else 1)
