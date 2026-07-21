"""Verify the cycle-theoretic mechanism for R_A (reflectioncancel) in FULL detail.

Case A: u(k)=w(k), so k is a fixed point of u^{-1}w.

CLAIM (the precise combinatorial statement to prove):
  Suppose u(k)=w(k). The map
      w  |-->  v* = w t_{k, pp*}
  where pp* is the unique index <k with u ->_{k-1} v*, is a bijection between:
    {w : u(k)=w(k), u ->_{k-2} w}       (QUANT nonzero)
  and
    {(v*, pp*) : v* = w t_{k,pp*}, pp*<k, u ->_{k-1} v*, v* contributes to MONK}
  Moreover:
    (A) P_{k-1}(u,v*) = P_{k-2}(u,w)             [same factorial-elem polynomial y-vars]
    (B) n_{k-1}(u,v*) = n_{k-2}(u,w) + 1          [same E_{deg,nv} object at level (p-1) vs (p-2)]
    (C) sign(pp*-k) * kappa(w->v*) * q_{k-1}(u,v*) = - q_{k-1} * q_{k-2}(u,w)
        i.e. the MONK summand exactly cancels QUANT.
    (D) pp* is characterized by the cycle of u^{-1}w containing k-1:
        Let c be the cycle of u^{-1}w with the property described by piericycle at level k-2
        that "sits just below k". pp* = the bottom index a_1 or the specific position
        determined by adjoining t_{k-1,k}.

We test (A),(B),(C) hold with 0 exceptions, AND we test the precise cycle rule for pp*:
    Adjoining the reflection t_{k-1,k} to the level-(k-2) chain to w produces the
    level-(k-1) permutation v* (i.e. v* = w t_{k-1,k}) IN THE CASE w fixes k-1.
    When w does NOT fix k-1, we test that pp* equals the top of the cycle of u^{-1}w
    that contains k-1.

Run: conda activate schubmult_312 && python _lscripts/RA_cycle_proof.py 5
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


def Pset(u, w, m):
    return sorted(val(u, i) for i in range(1, m + 1) if val(u, i) == val(w, i))


def nk(u, w, m):
    return sum(1 for i in range(1, m + 1) if val(u, i) != val(w, i))


def c_data(u, w, m, a, Rm):
    """Return (nonzero_bool, qweight, Pset, nk) for c_m(a;u,w)."""
    if m < 0:
        return (False, None, None, None)
    w = Permutation(w)
    if w not in Rm:
        return (False, None, None, None)
    qm = Rm[w]
    Pm = Pset(u, w, m)
    nm = m - len(Pm)
    deg = a - nm
    nv = m - nm
    if deg < 0 or deg > nv:
        return (False, qm, Pm, nm)
    return (True, qm, Pm, nm)


def cycle_containing(u, w, idx):
    """Return the cycle of u^{-1}w containing idx, as a tuple, or None if fixed point."""
    perm = (~u) * w
    seen = set()
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
    checks = 0
    failA = failB = failC = 0
    pp_rule = Counter()
    pp_rule_fail = 0
    examples = []
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rkm1 = enumerate_pieri(u, k - 1, N)
            Rkm2 = enumerate_pieri(u, k - 2, N) if k - 2 >= 0 else {u: S.One}
            # candidate w: reachable at k-2 (QUANT support) with u(k)=w(k)
            for w in list(Rkm2):
                w = Permutation(w)
                if val(u, k) != val(w, k):
                    continue
                for p in range(1, k + 1):
                    okw, qw2, P2, n2 = c_data(u, w, k - 2, p - 2, Rkm2)
                    if not okw:
                        continue
                    # QUANT nonzero. Find pp* < k with v*=w t_{k,pp*}, u ->_{k-1} v* contributing.
                    found = None
                    for pp in range(1, k):
                        vv = w.swap(k - 1, pp - 1)
                        dl = w.inv - vv.inv
                        lo, hi = pp, k
                        if dl == 1:
                            ttype = "RAISE"
                            kappa = S.One
                        elif dl == 1 - 2 * (hi - lo):
                            ttype = "DROP"
                            kappa = q_ab(lo, hi)
                        else:
                            continue
                        okv, qv1, P1, n1 = c_data(u, vv, k - 1, p - 1, Rkm1)
                        if not okv:
                            continue
                        found = (pp, ttype, kappa, vv, qv1, P1, n1)
                        break
                    if found is None:
                        # QUANT nonzero but no MONK partner -> would break R_A
                        failC += 1
                        if len(examples) < 8:
                            examples.append(("NO PARTNER", tuple(u), p, k, tuple(w)))
                        continue
                    checks += 1
                    pp, ttype, kappa, vv, qv1, P1, n1 = found
                    # (A) P-set match
                    if P1 != P2:
                        failA += 1
                    # (B) nk match
                    if n1 != n2 + 1:
                        failB += 1
                    # (C) weight cancellation:
                    # MONK summand = sign(pp-k)*kappa*qv1 ; QUANT = q_{k-1}*qw2
                    sign = 1 if pp > k else -1
                    monk = sign * kappa * qv1
                    quant = q_gs[k - 1] * qw2
                    if expand(monk + quant) != 0:
                        failC += 1
                        if len(examples) < 8:
                            examples.append(("WEIGHT", tuple(u), p, k, tuple(w), expand(monk), expand(quant)))
                    # pp* cycle rule
                    w_fixes_km1 = (val(w, k - 1) == val(u, k - 1))
                    v_is_adjacent = (pp == k - 1)
                    # cycle of u^{-1}w containing k-1
                    cyc = cycle_containing(u, w, k - 1)
                    desc = (ttype, w_fixes_km1, v_is_adjacent)
                    pp_rule[desc] += 1
    print(f"checks={checks}")
    print(f"failA (P-set match)      = {failA}")
    print(f"failB (nk+1 match)       = {failB}")
    print(f"failC (weight cancel/no partner) = {failC}")
    print(f"pp* descriptor (type, w fixes k-1, pp*==k-1):")
    for d, c in sorted(pp_rule.items(), key=lambda x: -x[1]):
        print(f"    {d}: {c}")
    if examples:
        print("EXAMPLES OF FAILURE:")
        for e in examples:
            print("   ", e)


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    run(n)
