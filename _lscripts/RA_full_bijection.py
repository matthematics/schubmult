"""Verify the FULL bijection for R_A: both directions + at-most-one.

Case A: u(k)=w(k). Claims to nail down for a rigorous proof:

  (I)  MONK sum has AT MOST ONE nonzero term (already known, re-verify).
  (II) MONK nonzero  <=>  QUANT nonzero  (i.e. u ->_{k-2} w).
       -- forward: if some pp<k gives u ->_{k-1} v* then u ->_{k-2} w.
       -- backward: if u ->_{k-2} w then exactly one pp<k contributes (RA_cycle_proof).
  (III) NO pp>k contributes in Case A (leveldichotomy(ii) consequence).
  (IV) The contributing pp* is characterized as follows:
        Let c = cycle of u^{-1}w containing k-1 (if any).
        Sub-case w fixes k-1 (k-1 NOT in any nontrivial cycle):  pp* = k-1, v*=w t_{k-1,k}.
        Sub-case w moves k-1: pp* = bottom index a_1 of the cycle of u^{-1}w containing k-1,
                              OR more precisely the position determined by splitting.
       We just record which explicit local rule matches, to state it in the paper.

Run: conda activate schubmult_312 && python _lscripts/RA_full_bijection.py 5
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


def cycle_containing(u, w, idx):
    perm = (~u) * w
    start = idx
    cyc = [start]
    cur = perm[start - 1]
    while cur != start:
        cyc.append(cur)
        cur = perm[cur - 1]
    if len(cyc) == 1:
        return None
    return tuple(cyc)


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    n_atmostone_fail = 0
    n_equiv_fail = 0
    n_ppgtk_fail = 0
    ppstar_rule = Counter()
    total = 0
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
            for w in cands:
                w = Permutation(w)
                if val(u, k) != val(w, k):
                    continue
                for p in range(1, k + 1):
                    # gather all contributing MONK terms
                    contrib_lt = []  # pp<k
                    contrib_gt = []  # pp>k
                    for pp in range(1, N + 1):
                        if pp == k:
                            continue
                        vv = w.swap(k - 1, pp - 1)
                        dl = w.inv - vv.inv
                        lo, hi = min(pp, k), max(pp, k)
                        if dl == 1 or dl == 1 - 2 * (hi - lo):
                            if c_nonzero(u, vv, k - 1, p - 1, Rkm1):
                                if pp < k:
                                    contrib_lt.append(pp)
                                else:
                                    contrib_gt.append(pp)
                    quant_nonzero = c_nonzero(u, w, k - 2, p - 2, Rkm2)
                    monk_nonzero = bool(contrib_lt or contrib_gt)
                    total += 1
                    # (I) at most one total
                    if len(contrib_lt) + len(contrib_gt) > 1:
                        n_atmostone_fail += 1
                    # (III) none with pp>k
                    if contrib_gt:
                        n_ppgtk_fail += 1
                    # (II) equivalence MONK<=>QUANT
                    if monk_nonzero != quant_nonzero:
                        n_equiv_fail += 1
                    # (IV) rule for pp*
                    if contrib_lt:
                        pp = contrib_lt[0]
                        w_fixes_km1 = (val(w, k - 1) == val(u, k - 1))
                        cyc = cycle_containing(u, w, k - 1)
                        if w_fixes_km1:
                            rule = ("wfix", pp == k - 1)
                        else:
                            # is pp* the bottom (a_1, the second-largest? smallest?) of the cycle?
                            # cycle_containing returns [k-1, perm(k-1), ...]. bottom = min of cycle.
                            cyc_min = min(cyc)
                            cyc_wo = [c for c in cyc if c != k - 1]
                            rule = ("wmove", pp == cyc_min, pp == min(cyc_wo) if cyc_wo else None, pp in cyc)
                        ppstar_rule[rule] += 1
    print(f"total Case-A (u(k)=w(k)) instances checked = {total}")
    print(f"(I)   MONK has >1 term            : {n_atmostone_fail} failures")
    print(f"(II)  MONK nonzero != QUANT nonzero: {n_equiv_fail} failures")
    print(f"(III) some pp>k contributes       : {n_ppgtk_fail} failures")
    print(f"(IV)  pp* rule histogram:")
    for r, c in sorted(ppstar_rule.items(), key=lambda x: -x[1]):
        print(f"    {r}: {c}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    run(n)
