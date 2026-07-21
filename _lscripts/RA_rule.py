"""Reveal the combinatorial RULE for pp* and the P-set match in R_A (Case A).

QUANT!=0  <=>  u ->_{k-2} w  (and u(k)=w(k)).  There is a unique down-neighbor pp*<k
with v*=w t_{k,pp*} contributing.  Print detailed structure to derive the rule:
  - u, w, v*, pp*, transition type
  - values u(k-1),u(k),w(k-1),w(k),w(pp*),v*(k),v*(pp*)
  - P_{k-2}(u,w) vs P_{k-1}(u,v*)
  - which position "enters" P when going k-2 -> k-1 for v*.

Run: conda activate schubmult_312 && python _lscripts/RA_rule.py 4
"""

import sys
from collections import Counter

from schubmult import *  # noqa: F401,F403
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, prod
from schubmult.symbolic.poly.variables import GeneratingSet

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


def c_nonzero(u, w, m, a, Rm):
    if m < 0:
        return False
    w = Permutation(w)
    if w not in Rm:
        return False
    Pm = Pset(u, w, m)
    nm = m - len(Pm)
    deg = a - nm
    nv = m - nm
    return 0 <= deg <= nv


def run(n, max_show=25):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    shown = 0
    rule_hist = Counter()
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rkm1 = enumerate_pieri(u, k - 1, N)
            Rkm2 = enumerate_pieri(u, k - 2, N) if k - 2 >= 0 else {u: S.One}
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
                    if not c_nonzero(u, w, k - 2, p - 2, Rkm2):
                        continue
                    found = None
                    for pp in range(1, k):
                        vv = w.swap(k - 1, pp - 1)
                        dl = w.inv - vv.inv
                        lo, hi = pp, k
                        if dl == 1:
                            ttype = "RAISE"
                        elif dl == 1 - 2 * (hi - lo):
                            ttype = "DROP"
                        else:
                            continue
                        if c_nonzero(u, vv, k - 1, p - 1, Rkm1):
                            found = (pp, ttype, vv)
                            break
                    if found is None:
                        continue
                    pp, ttype, vv = found
                    # Characterize the rule for pp*:
                    # v* = w t_{k,pp*}: swaps values at positions pp* and k.
                    # In Case A, w(k)=u(k). v*(k)=w(pp*), v*(pp*)=w(k)=u(k).
                    # P_{k-1}(u,v*) should equal P_{k-2}(u,w). The "entering" position for v* at level k-1:
                    #   is it position k-1 or pp*?
                    Pkm2_w = Pset(u, w, k - 2)
                    Pkm1_v = Pset(u, vv, k - 1)
                    # position that v* fixes among 1..k-1 that w didn't count in 1..k-2:
                    entering = [i for i in range(1, k) if val(u, i) == val(vv, i) and (i > k - 2 or val(u, i) != val(w, i))]
                    # descriptor
                    desc = (
                        ttype,
                        pp == k - 1,
                        val(w, k - 1) == val(u, k - 1),   # is k-1 fixed in w?
                        val(vv, k - 1) == val(u, k - 1),  # is k-1 fixed in v*?
                    )
                    rule_hist[desc] += 1
                    if shown < max_show:
                        print(f"u={tuple(u)} p={p} k={k} w={tuple(w)} pp*={pp} {ttype} v*={tuple(vv)}")
                        print(f"    u(k-1)={val(u,k-1)} u(k)={val(u,k)} | w(k-1)={val(w,k-1)} w(k)={val(w,k)} w(pp*)={val(w,pp)} | v*(k-1)={val(vv,k-1)} v*(k)={val(vv,k)} v*(pp*)={val(vv,pp)}")
                        print(f"    P_{{k-2}}(u,w)={Pkm2_w}  P_{{k-1}}(u,v*)={Pkm1_v}  entering(v*,<=k-1)={entering}")
                        shown += 1
    print(f"\nrule descriptor (type, pp*==k-1, w fixes k-1, v* fixes k-1) histogram:")
    for desc, cnt in sorted(rule_hist.items(), key=lambda x: -x[1]):
        print(f"    {desc}: {cnt}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    run(n)
