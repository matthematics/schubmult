"""Characterize the image set {v*} INTRINSICALLY and prove b=k via cycle structure.

We want to show: the R_A partner v* of w (u(k)=w(k), u->_{k-2} w) is exactly
characterized as a permutation with:
   (i)   u ->_{k-1} v*
   (ii)  u NOT ->_{k-2} v*
   (iii) in the cycle decomposition of u^{-1}v*, the (unique, by leveldichotomy/onedescent)
         cycle c* containing an index > k-2 has its TOP equal to k, i.e. the largest
         element of c* is exactly k.
And pieriknotkn1(u,v*,k-1) top b = top of c* = k.

Test that {v* : contributing partners} == {v : u->_{k-1} v, u NOT->_{k-2} v,
        cycle-top of the >(k-2) cycle == k}  and that the pieriknotkn1 top is that top.

Also test the DIRECT cycle claim proving b=k:
   Since u ->_{k-1} v* but NOT ->_{k-2} v*, some cycle c* of u^{-1}v* contains an
   index > k-2 (namely, since exactly one index of c* exceeds k-2, and to violate
   level k-2 that index must exceed... actually equals the top b*). Claim b*=k.
   Reason candidate: v* = w t_{pp*,k}, w(k)=u(k). So v*(k)=w(pp*), v*(pp*)=w(k)=u(k).
   Thus position k is where v* differs from u by receiving w(pp*); the value u(k)=w(k)
   moved to position pp*<k. So in u^{-1}v*, index k maps: (u^{-1}v*)(k)=u^{-1}(w(pp*)).
   Need: k is the LARGEST index in its cycle. Test directly.

Run: conda activate schubmult_312 && python _lscripts/RA_bk_cycle.py 5
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


def cycles(perm, maxidx):
    seen = set()
    out = []
    for s in range(1, maxidx + 1):
        if s in seen or perm[s - 1] == s:
            continue
        cyc = [s]
        seen.add(s)
        cur = perm[s - 1]
        while cur != s:
            cyc.append(cur)
            seen.add(cur)
            cur = perm[cur - 1]
        out.append(tuple(cyc))
    return out


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    tot = 0
    fail_ktop = 0
    fail_vk = 0
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
                    # cycle of u^{-1} v* containing an index > k-2
                    perm = (~u) * vv
                    cyc_list = cycles(perm, N)
                    # find the cycle(s) with an element > k-2 that isn't <=k-2-only
                    big_cycles = [c for c in cyc_list if max(c) > k - 2]
                    # by leveldichotomy there should be exactly one cycle preventing level k-2
                    # its top (max) should be k
                    tops = [max(c) for c in big_cycles]
                    # The one relevant: the cycle whose top > k-2 and which contains k or k-1
                    rel = [c for c in big_cycles if (k in c or (k - 1) in c)]
                    if not rel:
                        fail_ktop += 1
                        continue
                    if any(max(c) != k for c in rel):
                        fail_ktop += 1
                    # also v*(pp*) == u(k) and v* differs from u at position k
                    if val(vv, pp) != val(u, k):
                        fail_vk += 1
    print(f"contributing instances = {tot}")
    print(f"relevant cycle top != k        : {fail_ktop} failures")
    print(f"v*(pp*) != u(k)                : {fail_vk} failures")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    run(n)
