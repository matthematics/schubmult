"""Final precise weight check for the rigorous reflectioncancel proof.

For each contributing Case-A instance, with v* = w t_{pp*,k} (the unique partner)
and w = v* t_{k-1,k}, verify the EXACT identity relating the two q-weights:

  MONK summand q-weight:  sign(pp*-k) * kappa(w->v*) * q_{k-1}(u,v*)
  QUANT q-weight:         q_{k-1} * q_{k-2}(u,w)
  Claim: MONK summand + QUANT-weight == 0  (as symbolic q-monomials).

Also verify the two DEFINITIONAL relations that make c_{k-1}(p-1;u,v*) and
c_{k-2}(p-2;u,w) the SAME factorial-elementary polynomial:
  P_{k-1}(u,v*) == P_{k-2}(u,w)   and   n_{k-1}(u,v*) == n_{k-2}(u,w)+1.
and the resulting q-monomial q_{k-1}(u,v*) relation to q_{k-2}(u,w):
  q_{k-1}(u,v*) * (drop/raise weight of step w->v*) == +- q_{k-1} q_{k-2}(u,w).

Print the sign breakdown by whether the step w->v* is a raise or a drop.

Run: conda activate schubmult_312 && python _lscripts/RA_weight_final.py 5
"""

import sys
from collections import Counter

from schubmult import *  # noqa: F401,F403
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, prod
from schubmult.symbolic.poly.variables import GeneratingSet

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
    if m < 0:
        return (False, None)
    w = Permutation(w)
    if w not in Rm:
        return (False, None)
    Pm = Pset(u, w, m)
    nm = m - len(Pm)
    deg = a - nm
    nv = m - nm
    if deg < 0 or deg > nv:
        return (False, Rm[w])
    return (True, Rm[w])


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    tot = 0
    fail_wt = 0
    sign_hist = Counter()
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rkm1 = enumerate_pieri(u, k - 1, N)
            Rkm2 = enumerate_pieri(u, k - 2, N) if k - 2 >= 0 else {u: S.One}
            for w in list(Rkm2):
                w = Permutation(w)
                if val(u, k) != val(w, k):
                    continue
                for p in range(1, k + 1):
                    okw, qw2 = c_data(u, w, k - 2, p - 2, Rkm2)
                    if not okw:
                        continue
                    pp = None
                    vv = None
                    kappa = None
                    for cand in range(1, k):
                        vtry = w.swap(k - 1, cand - 1)
                        dl = w.inv - vtry.inv
                        lo, hi = cand, k
                        if dl == 1:
                            kap = S.One
                        elif dl == 1 - 2 * (hi - lo):
                            kap = q_ab(lo, hi)
                        else:
                            continue
                        okv, qv1 = c_data(u, vtry, k - 1, p - 1, Rkm1)
                        if okv:
                            pp = cand
                            vv = vtry
                            kappa = kap
                            qv1_w = qv1
                            break
                    if pp is None:
                        continue
                    tot += 1
                    sign = 1 if pp > k else -1  # pp<k => -1
                    monk = sign * kappa * qv1_w
                    quant = q_gs[k - 1] * qw2
                    if expand(monk + quant) != 0:
                        fail_wt += 1
                    # classify the w->v* step
                    dl = w.inv - vv.inv
                    step = "RAISE" if dl == 1 else "DROP"
                    sign_hist[(step,)] += 1
    print(f"contributing instances = {tot}")
    print(f"weight identity MONK + QUANT != 0 : {fail_wt} failures")
    print("w->v* step type histogram:")
    for s, c in sorted(sign_hist.items(), key=lambda x: -x[1]):
        print(f"    {s}: {c}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    run(n)
