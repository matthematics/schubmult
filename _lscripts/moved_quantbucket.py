"""For QUANT buckets in the moved case, verify the invariant/weight matching that the
quantum-term cancellation needs (so we can prove it WITHOUT citing reflectioncancel,
whose proof uses k fixed). Also record whether B==k in quant buckets and a* structure.

QUANT bucket: weight == q_{k-1}q_{k-2}(u,w), u->_{k-2}w. Single down t_{a*k}.
Check (n=5 AND n=6):
  (Q0) B == k  (C top is k)
  (Q1) P_{k-1}(u, v_dn) == P_{k-2}(u,w)      [set]
  (Q2) n_{k-1}(u, v_dn) == n_{k-2}(u,w)+1     via n_m=m-|P_m|
  (Q3) kappa_{a*k} q_{k-1}(u, v_dn) == q_{k-1}q_{k-2}(u,w)  [expand]
  (Q4) down summand kappa_{a*k} c_{k-1}(p-1;u,v_dn) == -(quantum) for a probe p  [expand, using E]
  (Q5) a* in C ; is a* a peeled neighbor (C-split leaves k a fixed pt or top-k piece)?
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
def Pset(u, w, m): return frozenset(val(u, i) for i in range(1, m+1) if val(u, i) == val(w, i))
def cyc_through(perm, x, N):
    c = [x]; cur = perm[x-1]
    while cur != x: c.append(cur); cur = perm[cur-1]
    return c
def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    tot = 0; f0 = f1 = f2 = f3 = 0; ain = 0
    for u in perms:
        u = Permutation(u)
        for k in range(3, n):  # need k-2>=1
            Rkm1 = enumerate_pieri(u, k-1, N); Rk = enumerate_pieri(u, k, N)
            Rkm2 = enumerate_pieri(u, k-2, N)
            cands = set(Rkm1) | set(Rk)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cands.add(v0.swap(k-1, pp-1))
            for w in cands:
                w = Permutation(w)
                if val(u, k) == val(w, k) or w in Rk: continue
                if w not in Rkm2: continue  # need quantum term nonzero
                quantwt = expand(q_ab(k-1, k)*Rkm2[w])
                downs = []
                for a in range(1, k):
                    vv = w.swap(a-1, k-1)
                    if (w.inv - vv.inv) in (1, 1-2*(k-a)) and vv in Rkm1:
                        kap = S.One if (w.inv - vv.inv) == 1 else q_ab(a, k)
                        if expand(kap*Rkm1[vv]) == quantwt:
                            downs.append((a, kap, vv))
                if len(downs) != 1: continue
                tot += 1
                a, kap, vv = downs[0]
                pi = (~u)*w; C = cyc_through(pi, k, N); B = max(C)
                if B != k: f0 += 1
                if a in C: ain += 1
                if Pset(u, vv, k-1) != Pset(u, w, k-2): f1 += 1
                nkm1 = (k-1) - len(Pset(u, vv, k-1)); nkm2 = (k-2) - len(Pset(u, w, k-2))
                if nkm1 != nkm2 + 1: f2 += 1
                if expand(kap*Rkm1[vv]) != quantwt: f3 += 1
    print(f"n={n}: quant buckets w/ unique down = {tot}")
    print(f"  (Q0) B==k                       : fails={f0}")
    print(f"  (Q1) P_{{k-1}}(u,v)==P_{{k-2}}(u,w) : fails={f1}")
    print(f"  (Q2) n_{{k-1}}==n_{{k-2}}+1        : fails={f2}")
    print(f"  (Q3) weight==q_{{k-1}}q_{{k-2}}     : fails={f3}")
    print(f"  a* in C                         : {ain}/{tot}")
if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
