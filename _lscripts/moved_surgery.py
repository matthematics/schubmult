"""Nail the cycle surgery between up t_{k,B} and down t_{a*,k}.
For each (u,w,k) moved case with a contributing up v_up=wt_{kB}:
  - print C (pi-cycle through k), B=max C, its indices >= k
  - a* (the down index), whether a* in C
  - structure of sigma_up=u^{-1}v_up cycles vs sigma_dn=u^{-1}v_dn
  - relation: pi = sigma_up t_{kB} = sigma_dn t_{a*k}
Report aggregate invariants:
  (S1) |C indices >= k| == 2  (exactly k and B)
  (S2) a* in C
  (S3) in sigma_up, k and B are in DIFFERENT cycles (split); in sigma_dn too
Print first 25 detailed cases.
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
def cycles_of(perm, N):
    seen = set(); cyc = []
    for i in range(1, N+1):
        if i in seen or perm[i-1] == i: continue
        c = [i]; seen.add(i); cur = perm[i-1]
        while cur != i: c.append(cur); seen.add(cur); cur = perm[cur-1]
        cyc.append(c)
    return cyc
def cycle_through(perm, x, N):
    c = [x]; cur = perm[x-1]
    while cur != x: c.append(cur); cur = perm[cur-1]
    return c
def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    s1 = s2 = s3 = 0; tot = 0; shown = 0
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
                        ups.append(b)
                if not ups: continue
                downs = []
                for a in range(1, k):
                    vv = w.swap(a-1, k-1)
                    if (w.inv - vv.inv) in (1, 1-2*(k-a)) and vv in Rkm1:
                        downs.append(a)
                if len(ups) != 1 or len(downs) != 1: 
                    tot += 1; continue
                tot += 1
                b = ups[0]; a = downs[0]
                C = cycle_through(pi, k, N); B = max(C)
                Cge = sorted([x for x in C if x >= k])
                if len(Cge) == 2 and Cge == [k, B]: s1 += 1
                if a in C: s2 += 1
                v_up = w.swap(k-1, b-1); v_dn = w.swap(a-1, k-1)
                su = (~u)*v_up; sd = (~u)*v_dn
                cu_k = cycle_through(su, k, N); cu_B = cycle_through(su, B, N)
                sd_k = cycle_through(sd, k, N)
                split_up = (B not in cu_k)
                if split_up: s3 += 1
                if shown < 25:
                    shown += 1
                    print(f"u={list(u)[:n]} w={list(w)[:n]} k={k} | C={C} B={B} Cge={Cge} | up b={b} dn a={a} a_in_C={a in C}")
                    print(f"    pi cycles: {cycles_of(pi,N)}")
                    print(f"    sigma_up cycles: {cycles_of(su,N)}  (k-cyc={cu_k}, B-cyc={cu_B})")
                    print(f"    sigma_dn cycles: {cycles_of(sd,N)}  (k-cyc={sd_k})")
    print(f"\nn={n}: cases(unique up&down)={s1 and tot}  tot_up_exists={tot}")
    print(f"  (S1) |C_>=k|==2 (k,B)      : {s1}")
    print(f"  (S2) a* in C               : {s2}")
    print(f"  (S3) sigma_up splits k|B   : {s3}")
if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
