"""Verify: for v* (the R_A partner), the UNIQUE b>k-1 from pieriknotkn1 at level k-1
is ALWAYS exactly b=k. This makes R_A's map w -> v* the inverse of pieriknotkn1's
map v* -> v* t_{k-1,b}, hence a bijection, and R_A follows from the ESTABLISHED lemma.

pieriknotkn1(u, v*, level k-1): since u ->_{k-1} v* and u NOT ->_{k-2} v*,
there is a unique b>k-1 with v* t_{k-1,b} <| v* and u ->_{k-2} v* t_{k-1,b}.
We enumerate all such b and check the set is exactly {k}.

Run: conda activate schubmult_312 && python _lscripts/RA_unique_b_is_k.py 5
"""

import sys

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


def reachable(u, w, k, N):
    return Permutation(w) in enumerate_pieri(u, k, N)


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


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    tot = 0
    fail_unique = 0
    fail_isk = 0
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
                    pp = None
                    vv = None
                    for cand in range(1, k):
                        vtry = w.swap(k - 1, cand - 1)
                        dl = w.inv - vtry.inv
                        lo, hi = cand, k
                        if dl == 1 or dl == 1 - 2 * (hi - lo):
                            if c_nonzero(u, vtry, k - 1, p - 1, Rkm1):
                                pp = cand
                                vv = vtry
                                break
                    if pp is None:
                        continue
                    tot += 1
                    # pieriknotkn1 at level k-1 for v*: enumerate all b>k-1
                    bs = []
                    for b in range(k, N + 1):
                        back = vv.swap(k - 2, b - 1)  # t_{k-1,b}
                        dl2 = vv.inv - back.inv
                        is_lup = (dl2 == 1) or (dl2 == 1 - 2 * (b - (k - 1)))
                        if is_lup and reachable(u, back, k - 2, N):
                            bs.append(b)
                    if len(bs) != 1:
                        fail_unique += 1
                    elif bs[0] != k:
                        fail_isk += 1
    print(f"contributing instances = {tot}")
    print(f"pieriknotkn1 gives non-unique b>k-1 : {fail_unique} failures")
    print(f"the unique b != k                   : {fail_isk} failures")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    run(n)
