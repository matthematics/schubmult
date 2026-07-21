"""Decisive check of the PROOF's algebra, independent of the unproven Lemmas 3.9-3.11.

The inductive proof claims that, extracting the coefficient of S^q_w from
    S_u E^q_{p,k} = (x_k - z_{k-p+1}) S_u E^q_{p-1,k-1} + S_u E^q_{p,k-1} + q_{k-1} S_u E^q_{p-2,k-2}
via the single-variable quantum Monk (Lemma 3.2) at i=k, gives exactly c_k(p;u,w), where
    c_m(a;u,w) = q_m(u,w) E_{a-n_m(u,w), m-n_m(u,w)}(y_{P_m(u,w)}; z)   (0 if not u ->_m w).

This script implements c_m FROM SCRATCH (Definition 1.1 enumerated by BFS + factorial elementary
in y_{P}, z) and rebuilds the coefficient
    coeff(u,p,k,w) = (y_{w(k)} - z_{k-p+1}) c_{k-1}(p-1;u,w)                      [diagonal]
        + sum_{pp != k}  sign(pp-k) * [q_{k,pp} if drop] * c_{k-1}(p-1; u, w t_{k,pp})  [Monk off-diag]
        + c_{k-1}(p;u,w)                                                          [E_{p,k-1}]
        + q_{k-1} c_{k-2}(p-2;u,w)                                                [quantum recursion]
and checks coeff(u,p,k,w) == c_k(p;u,w) as a symbolic identity for all u,p,k,w.

If this PASSES universally, the proof's algebra/bookkeeping is CORRECT and only the *justifying*
lemmas (why the terms organize as claimed) are wrong -- i.e. the theorem's proof needs correct
lemmas, not new algebra.  If it FAILS, the drafted proof has an actual algebraic error.

It also classifies, per (u,p,k,w), the four cases and tallies which Monk off-diagonal terms are
RAISE vs DROP and up (pp>k) vs down (pp<k), to reveal the true cancellation structure.

Run: conda activate schubmult_312 && python _lscripts/stress_test_proof_algebra.py 4
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
    """{w: q_k(u,w)} for all w with u ->_k w (Definition 1.1)."""
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
    """c_m(a;u,w) = q_m E_{a-n_m, m-n_m}(y_{P_m}; z), 0 unless u->_m w and n_m<=a<=m."""
    if m < 0:
        # E^q_{p,k}=0 for k<0 except the empty conventions handled by callers
        return S.Zero
    w = Permutation(w)
    if w not in Rm:
        return S.Zero
    qm = Rm[w]
    Pm = sorted(val(u, i) for i in range(1, m + 1) if val(u, i) == val(w, i))
    nm = m - len(Pm)
    deg = a - nm
    nv = m - nm  # = len(Pm)
    if deg < 0 or deg > nv:
        return S.Zero
    if nv == 0:
        return qm if deg == 0 else S.Zero
    yvars = [y[v] for v in Pm]
    ncoeff = nv + 1 - deg
    zvars = [z[i] for i in range(1, ncoeff + 1)]
    return qm * expand_func(FactorialElemSym(deg, nv, yvars, zvars))


def rebuild_coeff(u, p, k, w, N, Rk, Rkm1, Rkm2):
    u = Permutation(u)
    w = Permutation(w)
    j = k - p + 1
    total = S.Zero
    # diagonal
    total += (y[val(w, k)] - z[j]) * c_term(u, w, k - 1, p - 1, Rkm1)
    # Monk off-diagonal: sources v = w t_{k,pp}
    for pp in range(1, N + 1):
        if pp == k:
            continue
        v = w.swap(k - 1, pp - 1)
        dl = w.inv - v.inv  # length change v -> w
        sign = 1 if pp > k else -1
        lo, hi = min(k, pp), max(k, pp)
        if dl == 1:
            total += sign * c_term(u, v, k - 1, p - 1, Rkm1)
        elif dl == 1 - 2 * (hi - lo):
            total += sign * q_ab(lo, hi) * c_term(u, v, k - 1, p - 1, Rkm1)
    # E_{p,k-1} recursion term
    total += c_term(u, w, k - 1, p, Rkm1)
    # quantum recursion term
    total += q_gs[k - 1] * c_term(u, w, k - 2, p - 2, Rkm2)
    return total


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    total_checks = 0
    fails = []
    case_counts = defaultdict(int)
    for u in perms:
        for k in range(2, n):
            Rk = enumerate_pieri(u, k, N)
            Rkm1 = enumerate_pieri(u, k - 1, N)
            Rkm2 = enumerate_pieri(u, k - 2, N) if k - 2 >= 0 else {Permutation(u): S.One}
            # candidate w: union of supports and Monk-reachable
            cands = set(Rk) | set(Rkm1) | set(Rkm2)
            for v in list(Rkm1):
                for pp in range(1, N + 1):
                    if pp != k:
                        cands.add(v.swap(k - 1, pp - 1))
            for p in range(1, k + 1):
                for w in cands:
                    lhs = rebuild_coeff(u, p, k, w, N, Rk, Rkm1, Rkm2)
                    rhs = c_term(u, w, k, p, Rk)
                    total_checks += 1
                    if expand(lhs - rhs) != 0:
                        fails.append((tuple(u), p, k, tuple(w), str(expand(lhs - rhs))))
                    in_k = Permutation(w) in Rk
                    in_km1 = Permutation(w) in Rkm1
                    case_counts[(in_k, in_km1)] += 1
    print(f"n={n}: checked {total_checks} (u,p,k,w) coefficient identities; failures={len(fails)}")
    print(f"  case distribution (in_k, in_km1): {dict(case_counts)}")
    for f in fails[:20]:
        print("  ALGEBRA FAIL:", f)
    return len(fails) == 0


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    ok = run(n)
    print("RESULT:", "proof algebra CONSISTENT" if ok else "proof algebra HAS ERRORS")
    sys.exit(0 if ok else 1)
