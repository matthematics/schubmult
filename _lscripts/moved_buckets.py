"""Classify every Case-(ii) contributing Monk summand by q-bucket and confirm taxonomy.
For each moved case (u not->_k w, u(k)!=w(k)), group up/down summands by q-weight (expand key).
Report bucket types:
  DIAG: weight == q_{k-1}(u,w) expand  (Case i; may contain diagonal/straight) -> skip from (ii)
  QUANT: weight == q_{k-1}*q_{k-2}(u,w) AND u->_{k-2}w
  OTHER: everything else
For QUANT buckets: check (#up==0, #down==1).
For OTHER buckets: check (#up==1, #down==1) and up at b=B.
Also global: does sum over each bucket of kappa*c (expand) == 0 including quantum term where relevant.
"""
import sys
from collections import defaultdict
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, prod
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import FactorialElemSym
q_gs = GeneratingSet("q"); y = GeneratingSet("y"); z = GeneratingSet("z")
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
def Pset(u, w, m): return sorted(val(u, i) for i in range(1, m+1) if val(u, i) == val(w, i))
def cyc_through(perm, x, N):
    c = [x]; cur = perm[x-1]
    while cur != x: c.append(cur); cur = perm[cur-1]
    return c
def cval(pfrac, m, u, v):
    # c_m(pfrac; u, v) as elem-sym; returns (qweight_expr, Efactor) or None
    R = pfrac
    return None
def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    fq = fo = fbucket = 0; nq = no = 0
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rkm1 = enumerate_pieri(u, k-1, N); Rk = enumerate_pieri(u, k, N)
            Rkm2 = enumerate_pieri(u, k-2, N) if k >= 2 else {}
            cands = set(Rkm1) | set(Rk)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cands.add(v0.swap(k-1, pp-1))
            for w in cands:
                w = Permutation(w)
                if val(u, k) == val(w, k) or w in Rk: continue
                # collect up/down summands with (weight-expr)
                ups = []; downs = []
                for b in range(k+1, N+1):
                    vv = w.swap(k-1, b-1)
                    if (w.inv - vv.inv) in (1, 1-2*(b-k)) and vv in Rkm1:
                        kap = S.One if (w.inv - vv.inv) == 1 else q_ab(k, b)
                        ups.append((b, expand(kap*Rkm1[vv])))
                for a in range(1, k):
                    vv = w.swap(a-1, k-1)
                    if (w.inv - vv.inv) in (1, 1-2*(k-a)) and vv in Rkm1:
                        kap = S.One if (w.inv - vv.inv) == 1 else q_ab(a, k)
                        downs.append((a, expand(kap*Rkm1[vv])))
                diagw = expand(Rkm1[w]) if w in Rkm1 else None
                quantw = expand(q_ab(k-1, k)*Rkm2[w]) if (w in Rkm2) else None
                buckets = defaultdict(lambda: {"up": [], "dn": []})
                for b, wt in ups: buckets[wt]["up"].append(b)
                for a, wt in downs: buckets[wt]["dn"].append(a)
                pi = (~u)*w; C = cyc_through(pi, k, N); B = max(C)
                for wt, d in buckets.items():
                    if diagw is not None and wt == diagw:
                        continue  # Case (i)
                    if quantw is not None and wt == quantw:
                        nq += 1
                        if not (len(d["up"]) == 0 and len(d["dn"]) == 1): fq += 1
                    else:
                        no += 1
                        if not (len(d["up"]) == 1 and len(d["dn"]) == 1 and d["up"][0] == B): fo += 1
    print(f"n={n}: QUANT buckets={nq} OTHER buckets={no}")
    print(f"  QUANT (#up==0,#dn==1)          : fails={fq}")
    print(f"  OTHER (#up==1@B,#dn==1)        : fails={fo}")
if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
