"""Consolidated 0-fail verification of the pure-reflection pairing lemma (expand-based).

For every moved case (u NOT->_k w, u(k)!=w(k)), fixed p, and each PURE reflection q-group
(exactly one UP, one DOWN by moved_perp), verify ALL of:
  (V1) up index s == B == max(C)  and  kappa_up == 1   (up reflection t_{kB} length-increasing)
  (V2) down index r == pi^{-1}(k)   (predecessor of k in its cycle)   and  r < k
  (V3) C = cycle of pi through k has EXACTLY two indices >= k, namely k and B (B>k)
  (V4) P_{k-1}(u, v_up) == P_{k-1}(u, v_dn)  and  n_{k-1} equal   (=> identical E)
  (V5) expand(kappa_up * q_{k-1}(u,v_up)) == expand(kappa_dn * q_{k-1}(u,v_dn))  (same q-weight)
  (V6) 3-cycle relation:  t_{rk} t_{kB} == 3-cycle (r,k,B)  i.e. v_dn = v_up (r k B)-conjugate;
       concretely v_up = w t_{kB}, v_dn = w t_{rk}, and w t_{rk} = w t_{kB} t_{kB} t_{rk}.
"""
import sys
from collections import defaultdict
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
                results.setdefault(nperm, nqw)
                stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results

def Pset(u, w, m): return frozenset(val(u, i) for i in range(1, m+1) if val(u, i) == val(w, i))
def npart(u, w, m): return m - len(Pset(u, w, m))
def orbit(perm, start, N):
    c = [start]; cur = perm[start-1]
    while cur != start and cur not in c: c.append(cur); cur = perm[cur-1]
    return c

def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    tot = 0
    f1 = f2 = f3 = f4 = f5 = f6 = 0
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
                    if w in Rkm1: kind[str(expand(Rkm1[w]))].add('ds')
                    if w in Rkm2: kind[str(expand(q_gs[k-1]*Rkm2[w]))].add('quant')
                    for pp in range(1, N+1):
                        if pp == k: continue
                        vv = w.swap(k-1, pp-1); dl = w.inv - vv.inv
                        lo, hi = min(pp, k), max(pp, k)
                        if dl == 1: extra = S.One
                        elif dl == 1-2*(hi-lo): extra = q_ab(lo, hi)
                        else: continue
                        if vv not in Rkm1: continue
                        nm = npart(u, vv, k-1); deg = (p-1)-nm; nv = (k-1)-nm
                        if not (0 <= deg <= nv): continue
                        groups[str(expand(extra*Rkm1[vv]))].append((pp, extra, vv))
                    for qk, mem in groups.items():
                        if 'ds' in kind[qk] or 'quant' in kind[qk]: continue
                        ups = [m for m in mem if m[0] > k]; downs = [m for m in mem if m[0] < k]
                        if not (len(ups) == 1 and len(downs) == 1): continue
                        tot += 1
                        s, kap_up, v_up = ups[0]; r, kap_dn, v_dn = downs[0]
                        C = orbit(pi, k, N); B = max(C)
                        pinv_k = (~pi)[k-1]
                        if not (s == B and kap_up == S.One): f1 += 1
                        if not (r == pinv_k and r < k): f2 += 1
                        ge = [x for x in C if x >= k]
                        if not (len(ge) == 2 and set(ge) == {k, B} and B > k): f3 += 1
                        if not (Pset(u, v_up, k-1) == Pset(u, v_dn, k-1)
                                and npart(u, v_up, k-1) == npart(u, v_dn, k-1)): f4 += 1
                        if expand(kap_up*Rkm1[v_up]) != expand(kap_dn*Rkm1[v_dn]): f5 += 1
                        # V6: v_up = w t_{kB}, v_dn = w t_{rk}
                        if not (v_up == w.swap(k-1, B-1) and v_dn == w.swap(r-1, k-1)): f6 += 1
    print(f"n={n} pure pairs checked = {tot}")
    print(f"  (V1) s==B & kappa_up==1              fails = {f1}")
    print(f"  (V2) r==pi^-1(k) & r<k               fails = {f2}")
    print(f"  (V3) C has exactly {{k,B}} >= k, B>k  fails = {f3}")
    print(f"  (V4) P_{{k-1}},n_{{k-1}} agree up==dn   fails = {f4}")
    print(f"  (V5) q-weights agree (expand)        fails = {f5}")
    print(f"  (V6) v_up=wt_kB, v_dn=wt_rk           fails = {f6}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
