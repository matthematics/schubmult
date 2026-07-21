"""Verify the CLEAN closed form: up = t_{k,B}, down = t_{k-1,k}. Correct convention. n=5 AND n=6.

Moved case u NOT->_k w, u(k)!=w(k), whenever a contributing up-neighbour exists:
 (F1) exactly one contributing up-neighbour, at b=B=max C.
 (F2) exactly one contributing down-neighbour, at a = k-1.
 (F3) P_{k-1}(u,wt_{kB}) == P_{k-1}(u,wt_{k-1,k}) == P_{k-1}(u,w); equal n_{k-1}.
 (F4) kappa_{kB} q_{k-1}(u,wt_{kB}) == kappa_{k-1,k} q_{k-1}(u,wt_{k-1,k})  [expand].
 (F5) symmetric: whenever a contributing DOWN exists, a contributing UP exists (so counts match).
Also report: does a contributing down ALWAYS have a=k-1 (even if we found it via down-first)?
"""
import sys
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
def Pset(u, w, m): return frozenset(val(u, i) for i in range(1, m+1) if val(u, i) == val(w, i))
def orbit_from(perm, start):
    c = [start]; cur = perm[start-1]
    while cur != start: c.append(cur); cur = perm[cur-1]
    return c
def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    up_cases = 0; dn_cases = 0
    fF1 = fF2 = fF3 = fF4 = fF5 = 0
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
                        ups.append((b, kap, vv))
                downs = []
                for a in range(1, k):
                    vv = w.swap(a-1, k-1)
                    if (w.inv - vv.inv) in (1, 1-2*(k-a)) and vv in Rkm1:
                        kap = S.One if (w.inv - vv.inv) == 1 else q_ab(a, k)
                        downs.append((a, kap, vv))
                B = max(orbit_from(pi, k))
                if ups:
                    up_cases += 1
                    if not (len(ups) == 1 and ups[0][0] == B): fF1 += 1
                    if not (len(downs) == 1 and downs[0][0] == k-1): fF2 += 1
                    if downs:
                        b, ku, vu = ups[0]; a, kd, vd = downs[0]
                        if not (Pset(u, vu, k-1) == Pset(u, vd, k-1) == Pset(u, w, k-1)): fF3 += 1
                        if expand(ku*Rkm1[vu]) != expand(kd*Rkm1[vd]): fF4 += 1
                if downs:
                    dn_cases += 1
                    if not ups: fF5 += 1
    print(f"n={n}: up-cases={up_cases} down-cases={dn_cases}")
    print(f"  (F1) unique up at b=B                 : fails = {fF1}")
    print(f"  (F2) unique down at a=k-1            : fails = {fF2}")
    print(f"  (F3) P_{{k-1}} up==down==w             : fails = {fF3}")
    print(f"  (F4) kappa*q weights match           : fails = {fF4}")
    print(f"  (F5) down exists => up exists        : fails = {fF5}")
if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
