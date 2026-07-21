"""Examine j>1 cases (pi(k) != B): full cycle-surgery anatomy to find the intrinsic down rule.
Correct convention throughout."""
import sys
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, prod
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
def cyc_all(perm, N):
    seen = set(); cycs = []
    for x in range(1, N+1):
        if x in seen or perm[x-1] == x: continue
        c = orbit_from(perm, x); seen |= set(c); cycs.append(c)
    return cycs
def descent_u(u, cyc):
    # unique a_i in cyc (non-top) with u(a_i)>u(a_{i-1}) per paper; compute descents crudely:
    # normal form (a_p,...,a_1,b): b=max. orbit from b: [b, a_p, a_{p-1},...,a_1]. descent = a_i with u(a_i)>u(prev in reduced word)
    b = max(cyc); orb = orbit_from_perm_cycle(cyc, b)
    return b, orb
def orbit_from_perm_cycle(cyc, b):
    # cyc is a list orbit; rotate so it starts at b
    i = cyc.index(b); return cyc[i:]+cyc[:i]
def run(n, want=20):
    N = 2*n; perms = list(Permutation.all_permutations(n)); shown = 0
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
                if pi[k-1] == max(orbit_from(pi, k)): continue  # skip j==1 (pi(k)=B)
                ups = [b for b in range(k+1, N+1)
                       if (w.inv - w.swap(k-1, b-1).inv) in (1, 1-2*(b-k)) and w.swap(k-1, b-1) in Rkm1]
                if not ups: continue
                downs = [a for a in range(1, k)
                         if (w.inv - w.swap(a-1, k-1).inv) in (1, 1-2*(k-a)) and w.swap(a-1, k-1) in Rkm1]
                C = orbit_from(pi, k); B = max(C)
                orb = orbit_from_perm_cycle(C, B)  # [B, a_p, ..., a_1]
                # positions: orb = [B=orb0, orb1=a_p, ..., orb_last=a_1]; k somewhere
                kpos = orb.index(k)
                # u-values along orb
                uvals = [val(u, x) for x in orb]
                a = downs[0] if downs else None
                apos = orb.index(a) if a in orb else None
                print(f"u={tuple(u)} w={tuple(w)} k={k} B={B}")
                print(f"   orb[B,a_p..a_1]={orb}  u(orb)={uvals}  k at pos {kpos}")
                print(f"   pi(k)={pi[k-1]}  down a={a} at orb pos {apos}  (pi^-1(k)={(~pi)[k-1]})")
                shown += 1
                if shown >= want: return
if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
