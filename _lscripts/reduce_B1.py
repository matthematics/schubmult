"""Characterize the B1 residual (Case B1: u(k)!=w(k), u->_k w AND u->_{k-1} w):
   c_k(p;u,w) - DIAG - REC == MONK + QUANT.
Break MONK into up/down and raise/drop; break by term counts; check the sharp structure
(is it also a small-term quantum cancellation?), to state a precise B1 lemma.

Also verifies the key structural facts used in the writeup:
   c_k(p;u,w)   = q_k(u,w)   E^{(p-1)}          (E^{(a)}=E_{a-n_{k-1},(k-1)-n_{k-1}}(y_P;z))
   c_{k-1}(a)   = q_{k-1}(u,w) E^{(a)}
and reports the q_k(u,w) vs q_{k-1}(u,w) relationship.

Run: conda activate schubmult_312 && python _lscripts/reduce_B1.py 5
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


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    monk_term_hist = Counter()
    monk_type_hist = Counter()
    qk_vs_qkm1 = Counter()
    checks = 0
    fails = 0
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
                    if val(u, k) == val(w, k):
                        continue  # B only
                    in_k = w in Rk
                    in_km1 = w in Rkm1
                    if not in_km1:
                        continue  # B1 = {u(k)!=w(k), w in R_{k-1}}; note this forces w NOT in R_k
                    j = k - p + 1
                    ck = c_term(u, w, k, p, Rk)
                    DIAG = (y[val(w, k)] - z[j]) * c_term(u, w, k - 1, p - 1, Rkm1)
                    REC = c_term(u, w, k - 1, p, Rkm1)
                    QUANT = q_gs[k - 1] * c_term(u, w, k - 2, p - 2, Rkm2)
                    monk_terms = []
                    for pp in range(1, N + 1):
                        if pp == k:
                            continue
                        vv = w.swap(k - 1, pp - 1)
                        dl = w.inv - vv.inv
                        sign = 1 if pp > k else -1
                        lo, hi = min(k, pp), max(k, pp)
                        if dl == 1:
                            term = sign * c_term(u, vv, k - 1, p - 1, Rkm1)
                            ty = ("up" if pp > k else "dn", "RAISE")
                        elif dl == 1 - 2 * (hi - lo):
                            term = sign * q_ab(lo, hi) * c_term(u, vv, k - 1, p - 1, Rkm1)
                            ty = ("up" if pp > k else "dn", "DROP")
                        else:
                            continue
                        if expand(term) != 0:
                            monk_terms.append((ty, term))
                    MONK = sum((t[1] for t in monk_terms), S.Zero)
                    checks += 1
                    if expand((ck - DIAG - REC) - (MONK + QUANT)) != 0:
                        fails += 1
                    monk_term_hist[len(monk_terms)] += 1
                    for ty, _ in monk_terms:
                        monk_type_hist[ty] += 1
                    # c_k should be 0 in B1
                    ck_is_zero = expand(ck) == 0
                    qk_vs_qkm1["ck_zero" if ck_is_zero else "ck_NONZERO"] += 1
    print(f"n={n}: B1 checks={checks} residual_fails={fails}")
    print(f"  MONK #terms histogram: {dict(sorted(monk_term_hist.items()))}")
    print(f"  MONK type histogram:   {dict(monk_type_hist)}")
    print(f"  c_k==0 in B1? in_k count: {dict(qk_vs_qkm1)}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    run(n)
