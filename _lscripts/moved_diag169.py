"""Diagnose the 169 H1 failures: split the sub-conditions and print examples."""
import sys
from collections import Counter
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, prod
from schubmult.symbolic.poly.variables import GeneratingSet
q_gs = GeneratingSet("q")
def q_ab(a, b): return prod([q_gs[s] for s in range(a, b)])
def val(perm, pos): return perm[pos-1]
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
def cyc_through(perm, x, N):
    c = [x]; cur = perm[x-1]
    while cur != x and cur not in c: c.append(cur); cur = perm[cur-1]
    return c
def run(n, want=15):
    N = 2*n; perms = list(Permutation.all_permutations(n)); shown = 0
    cnt = Counter()
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
                if val(u, k) == val(w, k): continue
                if w in Rk: continue
                pi = (~u)*w
                ups = []
                for b in range(k+1, N+1):
                    vv = w.swap(k-1, b-1); dl = w.inv - vv.inv
                    if dl == 1 or dl == 1-2*(b-k):
                        if vv in Rkm1: ups.append(b)
                if not ups: continue
                C = cyc_through(pi, k, N); B = max(C)
                geqk = sorted(x for x in C if x >= k)
                nup = len(ups); bB = (ups == [B]); twoidx = (geqk == sorted({k, B})); pik = (pi[k-1] == B)
                key = (nup > 1, not bB, not twoidx, not pik)
                cnt[key] += 1
                if key != (False, False, False, False) and shown < want:
                    shown += 1
                    print(f"u={tuple(u)} w={tuple(w)} k={k} ups={ups} C={C} B={B} geqk={geqk} pi(k)={pi[k-1]}")
    print(f"\n(key=nup>1, b!=B, C-not-two-{{k,B}}, pi(k)!=B) histogram:")
    for kk, vv in sorted(cnt.items()):
        print(f"  {kk}: {vv}")
if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
