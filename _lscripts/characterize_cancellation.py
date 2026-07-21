"""Characterize the TRUE cancellation mechanism replacing the false Lemmas 3.9-3.11.

Uses the c_m building blocks (Definition 1.1, implemented independently) to split the rebuilt
coefficient into a RAISE part and a DROP part:

  RAISE(u,p,k,w) = (y_{w(k)}-z_{k-p+1}) c_{k-1}(p-1;u,w)
                 + sum_{pp!=k, Monk RAISE (dl=+1)} sign(pp-k) c_{k-1}(p-1;u,w t_{k,pp})
                 + c_{k-1}(p;u,w)
  DROP(u,p,k,w)  = sum_{pp!=k, Monk DROP} sign(pp-k) q_{k,pp} c_{k-1}(p-1;u,w t_{k,pp})
                 + q_{k-1} c_{k-2}(p-2;u,w)

We test, per case (in_k = u ->_k w, in_km1 = u ->_{k-1} w), the refined hypotheses that would
give a CORRECT proof:

  H1  (Case in_k & in_km1):      DROP == 0   and   RAISE == c_k(p;u,w)
  H2  (Case in_k & not in_km1):  RAISE == 0 and DROP == c_k(p;u,w)     [pure quantum, or mixed?]
  H3  (Case not in_k & in_km1):  RAISE == 0   and   DROP == 0
  H4  (Case not in_k & not in_km1): RAISE == 0 and DROP == 0

We tally, per case, how many satisfy DROP==0, RAISE==0, RAISE==c_k, DROP==c_k, total==c_k, so we
learn the actual clean statement to prove.

Run: conda activate schubmult_312 && python _lscripts/characterize_cancellation.py 4
"""

import sys
from collections import defaultdict

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


def c_term(u, w, m, a, Rm):
    if m < 0:
        return S.Zero
    w = Permutation(w)
    if w not in Rm:
        return S.Zero
    qm = Rm[w]
    Pm = sorted(val(u, i) for i in range(1, m + 1) if val(u, i) == val(w, i))
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


def split_coeff(u, p, k, w, N, Rkm1, Rkm2):
    u = Permutation(u)
    w = Permutation(w)
    j = k - p + 1
    raise_part = (y[val(w, k)] - z[j]) * c_term(u, w, k - 1, p - 1, Rkm1)
    raise_part += c_term(u, w, k - 1, p, Rkm1)
    drop_part = q_gs[k - 1] * c_term(u, w, k - 2, p - 2, Rkm2)
    for pp in range(1, N + 1):
        if pp == k:
            continue
        v = w.swap(k - 1, pp - 1)
        dl = w.inv - v.inv
        sign = 1 if pp > k else -1
        lo, hi = min(k, pp), max(k, pp)
        if dl == 1:
            raise_part += sign * c_term(u, v, k - 1, p - 1, Rkm1)
        elif dl == 1 - 2 * (hi - lo):
            drop_part += sign * q_ab(lo, hi) * c_term(u, v, k - 1, p - 1, Rkm1)
    return raise_part, drop_part


def isz(e):
    return expand(e) == 0


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    tally = defaultdict(lambda: defaultdict(int))
    for u in perms:
        for k in range(2, n):
            Rk = enumerate_pieri(u, k, N)
            Rkm1 = enumerate_pieri(u, k - 1, N)
            Rkm2 = enumerate_pieri(u, k - 2, N) if k - 2 >= 0 else {Permutation(u): S.One}
            cands = set(Rk) | set(Rkm1) | set(Rkm2)
            for v in list(Rkm1):
                for pp in range(1, N + 1):
                    if pp != k:
                        cands.add(v.swap(k - 1, pp - 1))
            for p in range(1, k + 1):
                for w in cands:
                    rp, dp = split_coeff(u, p, k, w, N, Rkm1, Rkm2)
                    ck = c_term(u, w, k, p, Rk)
                    in_k = Permutation(w) in Rk
                    in_km1 = Permutation(w) in Rkm1
                    case = {(True, True): "C1", (True, False): "C2", (False, True): "C3", (False, False): "C4"}[(in_k, in_km1)]
                    t = tally[case]
                    t["total"] += 1
                    if isz(dp):
                        t["DROP==0"] += 1
                    if isz(rp):
                        t["RAISE==0"] += 1
                    if isz(rp - ck):
                        t["RAISE==ck"] += 1
                    if isz(dp - ck):
                        t["DROP==ck"] += 1
                    if isz(rp + dp - ck):
                        t["RAISE+DROP==ck"] += 1
    for case in ["C1", "C2", "C3", "C4"]:
        t = tally[case]
        print(f"{case}: {dict(t)}")
    return tally


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    run(n)
