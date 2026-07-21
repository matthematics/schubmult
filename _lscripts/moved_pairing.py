"""Discover the pure-reflection up<->down pairing structure for the cycle-surgery proof.

Moved case u NOT->_k w, u(k)!=w(k). pi=u^{-1}w. For each fixed p and each pure-reflection
q-weight group (exactly one UP r=s>k, one DOWN r<k by moved_perp), print:
  - pi cycle C through k written from its top B: [B=c_0,c_1,...], j=index of k, a0=pi(k) succ.
  - the UP index s and DOWN index r
  - kappa_up, kappa_dn (are they q-drops or 1?)
  - the 3-cycle relating wt_{sk} and wt_{rk}:  t_{rk} t_{ks} = (r,k,s)
  - positions of s and r within C (or outside C)
Goal: read off closed structural rule for (s,r) from C to drive a pieriknotkn1-style proof.
"""
import sys
from collections import defaultdict
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, expand_func, prod
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import FactorialElemSym

y = GeneratingSet("y"); z = GeneratingSet("z"); q_gs = GeneratingSet("q")
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
                results.setdefault(nperm, nqw)
                stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results

def Pset(u, w, m): return sorted(val(u, i) for i in range(1, m+1) if val(u, i) == val(w, i))
def npart(u, w, m): return m - len(Pset(u, w, m))

def has_arrow(u, w, m, R): return Permutation(w) in R

def orbit(perm, start, N):
    c = [start]; cur = perm[start-1]
    while cur != start and cur not in c: c.append(cur); cur = perm[cur-1]
    return c

def run(n, want=30):
    N = 2*n; perms = list(Permutation.all_permutations(n)); shown = 0
    rulehits = 0; total = 0
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rk = enumerate_pieri(u, k, N); Rkm1 = enumerate_pieri(u, k-1, N); Rkm2 = enumerate_pieri(u, k-2, N)
            cands = set(Rk) | set(Rkm1) | set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cands.add(v0.swap(k-1, pp-1))
            for w in cands:
                w = Permutation(w)
                if val(u, k) == val(w, k): continue
                if w in Rk: continue
                pi = (~u)*w
                for p in range(1, k+1):
                    groups = defaultdict(list); kind = defaultdict(set)
                    ci_qw = Rkm1.get(w)
                    if ci_qw is not None and npart(u, w, k-1) <= p-1 <= (k-1):
                        pass
                    # mark ds / quant weights
                    if w in Rkm1:
                        kind[str(expand(Rkm1[w]))].add('ds')
                    if w in Rkm2:
                        kind[str(expand(q_gs[k-1]*Rkm2[w]))].add('quant')
                    for pp in range(1, N+1):
                        if pp == k: continue
                        vv = w.swap(k-1, pp-1); dl = w.inv - vv.inv
                        lo, hi = min(pp, k), max(pp, k)
                        if dl == 1: extra = S.One
                        elif dl == 1-2*(hi-lo): extra = q_ab(lo, hi)
                        else: continue
                        if vv not in Rkm1: continue
                        # admissibility at this p: deg in [0,nv]
                        nm = npart(u, vv, k-1); deg = (p-1)-nm; nv = (k-1)-nm
                        if not (0 <= deg <= nv): continue
                        qk = str(expand(extra*Rkm1[vv]))
                        groups[qk].append((pp, extra, vv))
                    for qk, mem in groups.items():
                        if 'ds' in kind[qk] or 'quant' in kind[qk]: continue
                        ups = [m for m in mem if m[0] > k]; downs = [m for m in mem if m[0] < k]
                        if not (len(ups) == 1 and len(downs) == 1): continue
                        total += 1
                        s, kap_up, v_up = ups[0]; r, kap_dn, v_dn = downs[0]
                        C = orbit(pi, k, N); B = max(C)
                        orb = orbit(pi, B, N); j = orb.index(k); a0 = pi[k-1]
                        # hypothesize rule: s == B, and r == pi^{-1}(k) or pi(k)?
                        pinv_k = (~pi)[k-1]
                        rule = (s == B)
                        if rule: rulehits += 1
                        if shown < want:
                            shown += 1
                            print(f"u={tuple(u)} w={tuple(w)} k={k} p={p} | s(up)={s} r(dn)={r} "
                                  f"B={B} a0=pi(k)={a0} pi^-1(k)={pinv_k} orbB={orb} j={j} "
                                  f"kap_up={'1' if kap_up==S.One else expand(kap_up)} "
                                  f"kap_dn={'1' if kap_dn==S.One else expand(kap_dn)} sEqB={s==B}")
    print(f"\ntotal pure pairs={total}  s==B holds in {rulehits}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
