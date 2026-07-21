"""Debug the 1629 (A)-failures: is w really in R_{k-1} or R_{k-2}, and what is C?"""
import sys
from collections import Counter
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, prod
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
def cyc_all(perm, N):
    seen = set(); cycs = []
    for x in range(1, N+1):
        if x in seen or perm[x-1] == x: continue
        c = [x]; seen.add(x); cur = perm[x-1]
        while cur != x: c.append(cur); seen.add(cur); cur = perm[cur-1]
        cycs.append(c)
    return cycs
def run(n, want=12):
    N = 2*n; perms = list(Permutation.all_permutations(n)); shown = 0; cnt = Counter()
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rk = enumerate_pieri(u, k, N); Rkm1 = enumerate_pieri(u, k-1, N); Rkm2 = enumerate_pieri(u, k-2, N)
            cands = set(Rkm1) | set(Rk) | set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cands.add(v0.swap(k-1, pp-1))
            for w in cands:
                w = Permutation(w)
                if val(u, k) == val(w, k): continue
                if w in Rk: continue
                inkm1 = w in Rkm1; inkm2 = w in Rkm2
                if not (inkm1 or inkm2): continue
                cnt[(inkm1, inkm2)] += 1
                if shown < want:
                    shown += 1
                    pi = (~u)*w
                    print(f"u={tuple(u)} w={tuple(w)} k={k} inR(k-1)={inkm1} inR(k-2)={inkm2} "
                          f"pi-cycles={cyc_all(pi, N)}")
    print(f"\n(inR(k-1), inR(k-2)) histogram among (A)-fails: {dict(cnt)}")
if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
