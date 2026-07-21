"""Pin down the EXACT cycle rule for pp* in the 'w moves k-1' subcase of R_A.

We know pp* is always IN the cycle c of u^{-1}w containing k-1.
Write c in piericycle normal form c = (a_p, ..., a_1, b) with b>k-2 the top,
a_i <= k-2, unique descent (Lemma onedescent).

We want: express pp* in terms of the cycle structure, so we can prove the
splitting u ->_{k-2} w  --adjoin t_{k-1,k}-->  u ->_{k-1} v* with v*=w t_{pp*,k}.

For each 'wmove' instance, print:
   cycle c (containing k-1), its normal form, position of k-1 in c, descent, pp*,
   and the relation v* = w t_{pp*,k}.

Also test a SHARPER claim:
   In ALL Case-A contributing instances, v* = w t_{pp*,k} and simultaneously
   the level-(k-2) chain to w, followed by the transposition t_{k-1,k}, equals
   a level-(k-1) chain to v*.  Concretely test:  u ->_{k-2} w  AND  w = v* t_{pp*,k}
   AND  there is a length-1-or-drop step  v* --t_{k-1,k}--> (something) ... 
   Simplest testable surrogate:  v* t_{k-1,k} relationship to w.

Run: conda activate schubmult_312 && python _lscripts/RA_ppstar_rule.py 5
"""

import sys
from collections import Counter

from schubmult import *  # noqa: F401,F403
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, prod
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


def run(n, show=30):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    shown = 0
    # Test claim: pp* = the position j in cycle c (containing k-1) such that
    #   w(pp*) is the value that, after swapping with position k, restores level k-1 chain.
    # Sharper testable surrogate claims:
    claim_top = Counter()   # pp* == b (top of cycle, >k-2)?  can't be since pp*<k... skip
    claim_vstar_wtkm1k = 0  # does v* == w t_{k-1,k} even when w moves k-1?
    tot = 0
    rule_match = Counter()
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
                    if not c_nonzero(u, w, k - 2, p - 2, Rkm2):
                        continue
                    # find pp*
                    pp = None
                    for cand in range(1, k):
                        vv = w.swap(k - 1, cand - 1)
                        dl = w.inv - vv.inv
                        lo, hi = cand, k
                        if dl == 1 or dl == 1 - 2 * (hi - lo):
                            if c_nonzero(u, vv, k - 1, p - 1, Rkm1):
                                pp = cand
                                break
                    if pp is None:
                        continue
                    tot += 1
                    w_fixes_km1 = (val(w, k - 1) == val(u, k - 1))
                    cyc = cycle_containing(u, w, k - 1)
                    if w_fixes_km1:
                        continue  # already handled: pp*=k-1
                    # w moves k-1: cyc is not None and contains k-1
                    # normal form: top = max(cyc)
                    top = max(cyc)
                    # rotate so top is last
                    ti = cyc.index(top)
                    nf = cyc[ti + 1:] + cyc[:ti + 1]  # (..., top)
                    # nf = (a_p,...,a_1,b) with b=top
                    body = nf[:-1]  # a_p..a_1  (these are <= k-2 in piericycle terms, but here at level k-2)
                    # position of k-1 in body
                    # Test: is pp* the element right "below" k-1 in the cycle order,
                    # i.e. w-image relation. Let's just test several candidate rules:
                    b = top
                    a1 = body[-1] if body else None   # a_1 (smallest-index end)
                    ap = body[0] if body else None    # a_p
                    km1_pos = body.index(k - 1) if (k - 1) in body else None
                    # candidate: pp* = a_1 (the bottom)
                    rules = {
                        "pp==a1": pp == a1,
                        "pp==ap": pp == ap,
                        "pp==min_cyc": pp == min(cyc),
                        "pp==w_of_something": False,
                    }
                    # Most robust: pp* is the position such that w(pp*) plays role.
                    # Record w(k-1), w(k), w(pp*), u(k-1), and which body elt = pp*
                    matched = tuple(sorted(r for r, ok in rules.items() if ok))
                    rule_match[matched] += 1
                    if shown < show:
                        print(f"u={tuple(u)} p={p} k={k} w={tuple(w)} pp*={pp}")
                        print(f"    cycle(k-1)={cyc} normal_form(a_p..a_1,b)={nf} body={body}")
                        print(f"    a_1={a1} a_p={ap} min_cyc={min(cyc)} | matched_rules={matched}")
                        shown += 1
    print(f"\nwmove instances = {tot - sum(1 for _ in [] )}")  # informational
    print("rule match histogram (which simple rules hold):")
    for r, c in sorted(rule_match.items(), key=lambda x: -x[1]):
        print(f"    {r}: {c}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    run(n)
