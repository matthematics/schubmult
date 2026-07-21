"""Test whether R_A reduces to the ALREADY-PROVEN Lemma pieriknotkn1.

Idea: In Case A (u(k)=w(k), u ->_{k-2} w, QUANT nonzero), the contributing MONK
neighbor is v* = w t_{pp*,k} with u ->_{k-1} v*.  Consider v* at LEVEL k:
   Is it true that  u ->_k v*  but  u NOT ->_{k-1} v*  is FALSE here...
Actually v* has u ->_{k-1} v* (contributes to c_{k-1}).  So instead consider the
mirror: apply pieriknotkn1 to the pair (u, v*) to relate level k-1 and k-2.

Precisely test the REVERSE-DIRECTION bijection that makes R_A a corollary of an
established uniqueness result:

  For (u,w) with u(k)=w(k) and u ->_{k-2} w:
    Claim R: there is a unique pp*<k with u ->_{k-1} (w t_{pp*,k})  [call it v*],
             AND this v* satisfies u ->_{k-1} v* but u NOT ->_{k-2} v*.
    Then pieriknotkn1 applied at level (k-1) to v* gives a UNIQUE b>k-1 with
         v* t_{k-1,b} <| v*  and  u ->_{k-2} v* t_{k-1,b}.
    Claim: that unique b is exactly k, and v* t_{k-1,k} = w.
  i.e.  the pair (pp* down-step from w) and (b=k up-step from v*) are inverse,
        so R_A's uniqueness = pieriknotkn1's uniqueness.

Test: (a) u ->_{k-1} v* and NOT u ->_{k-2} v*;
      (b) the unique b from pieriknotkn1(u, v*, level k-1) equals k;
      (c) v* t_{k-1,k} == w.

Run: conda activate schubmult_312 && python _lscripts/RA_via_pieriknotkn1.py 5
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
    failA = failB = failC = 0
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
                    # (a) u ->_{k-1} v* and NOT u ->_{k-2} v*
                    r_km1 = (vv in Rkm1)
                    r_km2 = reachable(u, vv, k - 2, N)
                    if not (r_km1 and not r_km2):
                        failA += 1
                    # (b) pieriknotkn1 at level k-1 for v*: unique b>k-1 with v* t_{k-1,b} <| v*,
                    #     u ->_{k-2} v* t_{k-1,b}. Test that b=k works.
                    #     v* t_{k-1,k}:
                    back = vv.swap(k - 2, k - 1)  # t_{k-1,k}
                    dl2 = vv.inv - back.inv
                    is_lup = (dl2 == 1) or (dl2 == 1 - 2 * (k - (k - 1)))
                    if not (is_lup and reachable(u, back, k - 2, N)):
                        failB += 1
                    # (c) v* t_{k-1,k} == w
                    if back != w:
                        failC += 1
    print(f"contributing instances = {tot}")
    print(f"(a) u ->_{{k-1}} v* & NOT ->_{{k-2}} v* : {failA} failures")
    print(f"(b) v* t_{{k-1,k}} <| v* & ->_{{k-2}}   : {failB} failures")
    print(f"(c) v* t_{{k-1,k}} == w                  : {failC} failures")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    run(n)
