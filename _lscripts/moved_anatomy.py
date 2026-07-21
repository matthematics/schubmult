"""Dump full anatomy of ONE moved case with an up-neighbour: cycles of pi, pi t_kB, pi t_rk,
and where the down-partner r lives (in C or a spectator?), plus spectator validity at k-1."""
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
def cyc_all(perm, N):
    seen = set(); cycs = []
    for x in range(1, N+1):
        if x in seen or perm[x-1] == x: continue
        c = [x]; seen.add(x); cur = perm[x-1]
        while cur != x: c.append(cur); seen.add(cur); cur = perm[cur-1]
        cycs.append(c)
    return cycs
def valid_level(cyc, m):
    big = [x for x in cyc if x > m]
    return len(big) == 1  # exactly one index > m (its top)
def run(n, maxshow=6):
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
                ups = [b for b in range(k+1, N+1)
                       if (w.swap(k-1, b-1).inv - w.inv in (1, 1-2*(b-k))) and w.swap(k-1, b-1) in Rkm1]
                if not ups: continue
                downs = [a for a in range(1, k)
                         if (w.swap(a-1, k-1).inv - w.inv in (1, 1-2*(k-a))) and w.swap(a-1, k-1) in Rkm1]
                allcyc = cyc_all(pi, N)
                C = next(c for c in allcyc if k in c); B = max(C)
                spec = [c for c in allcyc if k not in c]
                spec_ok_km1 = all(valid_level(c, k-1) for c in spec)
                b = ups[0]; r = downs[0] if downs else None
                pi_up = cyc_all((~u)*w.swap(k-1, b-1), N)
                pi_dn = cyc_all((~u)*w.swap(r-1, k-1), N) if r else None
                print(f"u={tuple(u)} w={tuple(w)} k={k} B={B}")
                print(f"   pi cycles = {allcyc}   C={C}  spectators={spec}")
                print(f"   spectators valid at k-1? {spec_ok_km1}")
                print(f"   up b={b}: pi t_kB cycles = {pi_up}")
                print(f"   downs={downs}  r={r}: pi t_rk cycles = {pi_dn}   r in C? {r in C if r else None}")
                shown += 1
                if shown >= maxshow: return
if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
