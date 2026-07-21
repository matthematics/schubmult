"""Identify the special down-neighbor pp* in the Case-A residual, to nail a bijective proof.

For each Case-A (u,p,k,w) with QUANT != 0, record:
  - pp* (the unique nonzero MONK down-neighbor position)
  - transition type of w -> v=w t_{k,pp*} (RAISE/DROP)
  - whether pp* == k-1
  - structural descriptor: position in u of the value w(k); the ->_{k-2} path info
Goal: a closed rule for pp* and a value identity
     [q if drop] c_{k-1}(p-1;u, w t_{k,pp*}) == q_{k-1} c_{k-2}(p-2;u,w).

Run: conda activate schubmult_312 && python _lscripts/find_ppstar.py 5
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
    ppstar_rule = Counter()     # pp* - (k-1) offset
    ppstar_desc = Counter()     # (type, pp*==k-1?)
    wk_vs_wppstar = Counter()   # relation of w(k) and w(pp*)
    total = 0
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
                    QUANT = q_gs[k - 1] * c_term(u, w, k - 2, p - 2, Rkm2)
                    if expand(QUANT) == 0:
                        continue
                    # find the unique nonzero down MONK neighbor
                    found = None
                    for pp in range(1, k):
                        vv = w.swap(k - 1, pp - 1)
                        dl = w.inv - vv.inv
                        lo, hi = pp, k
                        if dl == 1:
                            cval = c_term(u, vv, k - 1, p - 1, Rkm1)
                            ttype = "RAISE"
                        elif dl == 1 - 2 * (hi - lo):
                            cval = c_term(u, vv, k - 1, p - 1, Rkm1)
                            ttype = "DROP"
                        else:
                            continue
                        if expand(cval) != 0:
                            found = (pp, ttype, vv)
                            break
                    total += 1
                    if found is None:
                        ppstar_desc[("NONE_FOUND", None)] += 1
                        continue
                    pp, ttype, vv = found
                    ppstar_rule[pp - (k - 1)] += 1
                    ppstar_desc[(ttype, pp == k - 1)] += 1
                    # relation: is w(pp*) the value that "wraps"?  compare w(pp*),w(k),u(k-1)
                    rel = []
                    if val(w, pp) == val(u, k - 1):
                        rel.append("w(pp*)=u(k-1)")
                    if val(w, k) == val(u, k):
                        rel.append("w(k)=u(k)")
                    wk_vs_wppstar[tuple(rel)] += 1
    print(f"n={n}: nonzero-QUANT Case-A instances={total}")
    print(f"  pp* - (k-1) offset histogram: {dict(sorted(ppstar_rule.items()))}")
    print(f"  (type, pp*==k-1) histogram:   {dict(ppstar_desc)}")
    print(f"  value relations:              {dict(wk_vs_wppstar)}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    run(n)
