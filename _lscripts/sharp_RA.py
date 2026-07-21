"""Verify the SHARP Case-A claim discovered by dissection:
  (i)  at most ONE MONK neighbor is nonzero (in Case A, u(k)=w(k));
  (ii) MONK == -QUANT  (so R_A = MONK+QUANT = 0);
  (iii) identify the special neighbor pp* and its transition type.

Also probe Case B2 / B1 residual term-counts to gauge their tractability.

Run: conda activate schubmult_312 && python _lscripts/sharp_RA.py 5
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


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    fail_i = 0
    fail_ii = 0
    checks = 0
    monk_count_hist = Counter()
    special_type = Counter()
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rk = enumerate_pieri(u, k, N)
            Rkm1 = enumerate_pieri(u, k - 1, N)
            Rkm2 = enumerate_pieri(u, k - 2, N) if k - 2 >= 0 else {u: S.One}
            cands = set(Rk) | set(Rkm1) | set(Rkm2)
            for v in list(Rkm1):
                for pp in range(1, N + 1):
                    if pp != k:
                        cands.add(v.swap(k - 1, pp - 1))
            for p in range(1, k + 1):
                for w in cands:
                    w = Permutation(w)
                    if val(u, k) != val(w, k):
                        continue
                    checks += 1
                    monk_terms = []
                    for pp in range(1, N + 1):
                        if pp == k:
                            continue
                        vv = w.swap(k - 1, pp - 1)
                        dl = w.inv - vv.inv
                        sign = 1 if pp > k else -1
                        lo, hi = min(k, pp), max(k, pp)
                        if dl == 1:
                            cval = c_term(u, vv, k - 1, p - 1, Rkm1)
                            term = sign * cval
                            ttype = ("up" if pp > k else "dn", "RAISE")
                        elif dl == 1 - 2 * (hi - lo):
                            cval = c_term(u, vv, k - 1, p - 1, Rkm1)
                            term = sign * q_ab(lo, hi) * cval
                            ttype = ("up" if pp > k else "dn", "DROP")
                        else:
                            continue
                        if expand(term) != 0:
                            monk_terms.append((pp, ttype, term))
                    monk_count_hist[len(monk_terms)] += 1
                    if len(monk_terms) > 1:
                        fail_i += 1
                    if len(monk_terms) == 1:
                        special_type[monk_terms[0][1]] += 1
                    MONK = sum((t[2] for t in monk_terms), S.Zero)
                    QUANT = q_gs[k - 1] * c_term(u, w, k - 2, p - 2, Rkm2)
                    if expand(MONK + QUANT) != 0:
                        fail_ii += 1
    print(f"n={n}: Case-A checks={checks}")
    print(f"  (i)  #nonzero MONK terms histogram: {dict(sorted(monk_count_hist.items()))}")
    print(f"       >1 term violations: {fail_i}")
    print(f"  (ii) MONK+QUANT!=0 failures: {fail_ii}")
    print(f"  (iii) special neighbor type when exactly 1: {dict(special_type)}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    run(n)
