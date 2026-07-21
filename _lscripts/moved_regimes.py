"""Verify the regime separation needed for the writeup.

Moved case: u NOT->_k w, u(k)!=w(k). pi=u^{-1}w, C=cycle through k, B=max(C).
Claims to verify (0-fail expected at n=5):

(A) QUANT nonzero (u ->_{k-2} w)  =>  B=k.   [q_{k-1}q_{k-2}(u,w) weight present]
(B) B=k world: every contributing Monk neighbour wt_{rk} (u->_{k-1} wt_{rk}, wt_{rk}<|w,
    i.e. appears in the recursion) has q-weight in { q_{k-1}(u,w), q_{k-1}q_{k-2}(u,w) }.
    And NO contributing up-neighbour t_{kb} (b>k) with weight != q_{k-1}(u,w).
(C) B>k world: DIAG=STRAIGHT=QUANT=0 (u NOT->_{k-1}w and u NOT->_{k-2}w).
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
                results.setdefault(nperm, nqw)
                stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results

def orbit(perm, start, N):
    c = [start]; cur = perm[start-1]
    while cur != start: c.append(cur); cur = perm[cur-1]
    return c

def qweight_here(u, v, kk, N):
    R = enumerate_pieri(u, kk, N)
    return expand(R[v]) if v in R else None

def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    failA = failB = failC = 0; nBk = nBgt = 0
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rkm1 = enumerate_pieri(u, k-1, N)
            Rkm2 = enumerate_pieri(u, k-2, N) if k >= 2 else {}
            Rk = enumerate_pieri(u, k, N)
            cset = set(Rkm1)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cset.add(v0.swap(k-1, pp-1))
            for w in cset:
                w = Permutation(w)
                if val(u, k) == val(w, k): continue
                if w in Rk: continue
                pi = (~u)*w; C = orbit(pi, k, N); B = max(C)
                w_in_km1 = w in Rkm1     # DIAG/STRAIGHT nonzero
                w_in_km2 = w in Rkm2     # QUANT nonzero
                q_diag = expand(Rkm1[w]) if w_in_km1 else None
                q_quant = expand(q_gs[k-1]*q_gs[k-2]*Rkm2[w]) if w_in_km2 else None
                # (A)
                if w_in_km2 and B != k: failA += 1
                # (C)
                if B > k:
                    nBgt += 1
                    if w_in_km1 or w_in_km2: failC += 1
                if B == k:
                    nBk += 1
                    # (B): all contributing Monk neighbours weight in {q_diag, q_quant}
                    for r in range(1, N+1):
                        if r == k: continue
                        vv = w.swap(k-1, r-1); dl = w.inv - vv.inv
                        lo, hi = min(r, k), max(r, k)
                        if dl == 1: kap = S.One
                        elif dl == 1-2*(hi-lo): kap = q_ab(lo, hi)
                        else: continue
                        if vv not in Rkm1: continue
                        wt = expand(kap*Rkm1[vv])
                        okset = set()
                        if q_diag is not None: okset.add(q_diag)
                        if q_quant is not None: okset.add(q_quant)
                        if wt not in okset:
                            failB += 1
    print(f"n={n}  B=k world={nBk}  B>k world={nBgt}")
    print(f"  (A) QUANT nonzero but B!=k : fails = {failA}")
    print(f"  (B) B=k Monk weight not in {{q_diag,q_quant}} : fails = {failB}")
    print(f"  (C) B>k but DIAG/STRAIGHT/QUANT nonzero : fails = {failC}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
