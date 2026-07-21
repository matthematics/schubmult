"""Dissect the Case A residual R_A: MONK + QUANT == 0  (when u(k)=w(k)).

For each (u,p,k,w) with u(k)=w(k), print every nonzero contributing term:
  - MONK neighbors pp!=k: v=w t_{k,pp}, direction up/down, raise/drop, dropcoeff, c_{k-1}(p-1;u,v)
  - QUANT: q_{k-1} c_{k-2}(p-2;u,w)
and show they sum to zero, revealing the pairing (which term kills which).

Goal: find an explicit sign-reversing involution / pairing that proves R_A=0.

Run: conda activate schubmult_312 && python _lscripts/dissect_RA.py 4
"""

import sys

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


def dissect(n, max_show=12):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    shown = 0
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
                        continue  # Case A only
                    terms = []  # (label, symbolic value)
                    for pp in range(1, N + 1):
                        if pp == k:
                            continue
                        vv = w.swap(k - 1, pp - 1)
                        dl = w.inv - vv.inv
                        sign = 1 if pp > k else -1
                        lo, hi = min(k, pp), max(k, pp)
                        if dl == 1:
                            cval = c_term(u, vv, k - 1, p - 1, Rkm1)
                            if expand(cval) != 0:
                                terms.append((f"MONK pp={pp} {'up' if pp>k else 'dn'} RAISE", sign * cval, tuple(vv)))
                        elif dl == 1 - 2 * (hi - lo):
                            cval = c_term(u, vv, k - 1, p - 1, Rkm1)
                            if expand(cval) != 0:
                                terms.append((f"MONK pp={pp} {'up' if pp>k else 'dn'} DROP q[{lo},{hi})", sign * q_ab(lo, hi) * cval, tuple(vv)))
                    qval = q_gs[k - 1] * c_term(u, w, k - 2, p - 2, Rkm2)
                    if expand(qval) != 0:
                        terms.append(("QUANT", qval, tuple(w)))
                    if not terms:
                        continue
                    total = expand(sum(t[1] for t in terms))
                    if shown < max_show and len(terms) >= 1:
                        print(f"\nu={tuple(u)} p={p} k={k} w={tuple(w)}  (#terms={len(terms)}, sum={'0' if total==0 else 'NONZERO!'})")
                        for lab, valexpr, vt in terms:
                            print(f"    {lab:36s} v={vt}  val={expand(valexpr)}")
                        shown += 1
    print(f"\n(shown {shown} Case-A instances)")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    dissect(n)
