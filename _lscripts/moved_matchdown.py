"""When an up-neighbour exists (Case ii pure), characterize the MATCHING down.
Test several closed-form candidates for the matching down index a*:
  (C1) a* = k-1
  (C2) a* = min contributing down a
  (C3) a* = pi(B)      (image of B=maxC under pi=u^{-1}w)
  (C4) a* = the unique down whose q-weight (expand) equals the up's weight
Also: how many contributing downs are there when up exists? distribution.
"""
import sys
from collections import Counter
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, prod
from schubmult.symbolic.poly.variables import GeneratingSet
q_gs = GeneratingSet("q")
def q_ab(a, b): return prod([q_gs[s] for s in range(a, b)])
def val(p, i): return p[i-1]
def enumerate_pieri(u, k, N):
    u = Permutation(u); results = {u: S.One}
    stack = [(u, frozenset(), N+1, S.One, u.inv)]; seen = set()
    while stack:
        perm, used_a, last_b, qw, clen = stack.pop(); key = (perm, used_a, last_b)
        if key in seen: continue
        seen.add(key)
        for a in range(1, k+1):
            if a in used_a: continue
            for b in range(k+1, N+1):
                if b > last_b: continue
                nperm = perm.swap(a-1, b-1); d = nperm.inv - clen
                if d == 1: nqw = qw
                elif d == -2*(b-a)+1: nqw = qw*q_ab(a, b)
                else: continue
                results.setdefault(nperm, nqw); stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results
def orbit_from(perm, start):
    c = [start]; cur = perm[start-1]
    while cur != start: c.append(cur); cur = perm[cur-1]
    return c
def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    up_cases = 0
    c1 = c2 = c3 = 0
    ndowns_dist = Counter()
    match_a_minus_k = Counter()  # (a* - k) for the matching down
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rkm1 = enumerate_pieri(u, k-1, N); Rk = enumerate_pieri(u, k, N)
            cands = set(Rkm1) | set(Rk)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cands.add(v0.swap(k-1, pp-1))
            for w in cands:
                w = Permutation(w)
                if val(u, k) == val(w, k) or w in Rk: continue
                pi = (~u)*w
                ups = []
                for b in range(k+1, N+1):
                    vv = w.swap(k-1, b-1)
                    if (w.inv - vv.inv) in (1, 1-2*(b-k)) and vv in Rkm1:
                        kap = S.One if (w.inv - vv.inv) == 1 else q_ab(k, b)
                        ups.append((b, expand(kap*Rkm1[vv])))
                if not ups: continue
                downs = []
                for a in range(1, k):
                    vv = w.swap(a-1, k-1)
                    if (w.inv - vv.inv) in (1, 1-2*(k-a)) and vv in Rkm1:
                        kap = S.One if (w.inv - vv.inv) == 1 else q_ab(a, k)
                        downs.append((a, expand(kap*Rkm1[vv])))
                up_cases += 1
                ndowns_dist[len(downs)] += 1
                upwt = ups[0][1]
                # matching down = one whose weight equals up weight
                matching = [a for (a, wt) in downs if wt == upwt]
                if len(matching) == 1:
                    astar = matching[0]
                    match_a_minus_k[astar - k] += 1
                    if astar == k-1: c1 += 1
                    if downs and astar == min(a for a, _ in downs): c2 += 1
                    if astar == pi[max(orbit_from(pi, k))-1]: c3 += 1
    print(f"n={n}: up-cases={up_cases}")
    print(f"  #downs distribution: {dict(sorted(ndowns_dist.items()))}")
    print(f"  matching down (a*-k) distribution: {dict(sorted(match_a_minus_k.items()))}")
    print(f"  (C1) a*==k-1        : {c1}/{up_cases}")
    print(f"  (C2) a*==min down a : {c2}/{up_cases}")
    print(f"  (C3) a*==pi(B)      : {c3}/{up_cases}")
if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
