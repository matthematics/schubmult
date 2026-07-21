"""Intrinsic characterization of the unique up/down neighbours in a pure group, for the proof.

We do NOT assume s=B or r=pi^-1(k). Instead verify the INTRINSIC structural facts a
pieriknotkn1-style proof will use (all expand-checked, must pass n=5 AND n=6):

Fix moved case (u NOT->_k w, u(k)!=w(k)), pi=u^{-1}w, C=cycle thru k, B=max(C).
Contributing Monk neighbour: pp!=k with w'=wt_{pp,k}, w'<|w (raise or drop), u->_{k-1} w',
and admissible at p (0<=deg<=nv). PURE group = q-weight not equal to any diag/straight/quant wt.

Verify:
 (I1) EVERY contributing pp lies in C (pp in cycle of pi thru k).
 (I2) For a contributing pp, the reflection t_{pp,k} SPLITS C into two cycles, one containing
      k and one containing B (i.e. it separates the two >=k indices). Equivalently pp is on the
      k-side... test: in cycle C, both k and B present; t_{pp,k} applied to pi splits so k and B
      land in different cycles. Count how many pp achieve this split (should relate to pairing).
 (I3) The UNIQUE up index s (s>k) satisfies: s is the top B  OR ... just record what s is
      relative to C's orbit; and the UNIQUE down index r (<k). Record (orbit-pos(s), orbit-pos(r)).
 (I4) The map w' = w t_{s,k}  <->  w'' = w t_{r,k} is the 3-cycle exchange: w t_{r,k} =
      (w t_{s,k}) t_{s,k} t_{r,k}, and t_{s,k} t_{r,k} = 3-cycle (s,k,r). Record whether
      pi t_{s,k} and pi t_{r,k} give the SAME cycle type on C (both split C into {k-piece},{B-piece}).
"""
import sys
from collections import defaultdict, Counter
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
def cyc_through(perm, x, N):
    c = [x]; cur = perm[x-1]
    while cur != x and cur not in c: c.append(cur); cur = perm[cur-1]
    return set(c)

def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    tot = 0; fI1 = 0; fI2 = 0
    s_pos = Counter(); r_pos = Counter(); split_same = Counter()
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
                C = cyc_through(pi, k, N); B = max(C)
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
                        groups[str(expand(extra*Rkm1[vv]))].append(pp)
                    for qk, mem in groups.items():
                        if 'ds' in kind[qk] or 'quant' in kind[qk]: continue
                        ups = [pp for pp in mem if pp > k]; downs = [pp for pp in mem if pp < k]
                        if not (len(ups) == 1 and len(downs) == 1): continue
                        tot += 1
                        s = ups[0]; r = downs[0]
                        if s not in C: fI1 += 1
                        if r not in C: fI1 += 1
                        # I2: does t_{s,k} split C separating k and B?
                        pis = pi.swap(s-1, k-1); ck_s = cyc_through(pis, k, N)
                        split_s = (B not in ck_s)
                        pir = pi.swap(r-1, k-1); ck_r = cyc_through(pir, k, N)
                        split_r = (B not in ck_r)
                        if not (split_s and split_r): fI2 += 1
                        # positions in orbit from B
                        orb = [B]; cur = pi[B-1]
                        while cur != B and cur not in orb: orb.append(cur); cur = pi[cur-1]
                        s_pos[('s', 'isB' if s == B else 'other')] += 1
                        r_pos[('r', orb.index(r) if r in orb else -1)] += 1
    print(f"n={n} pure pairs = {tot}")
    print(f"  (I1) up or down index not in C : fails = {fI1}")
    print(f"  (I2) t_sk or t_rk fails to split C separating k,B : fails = {fI2}")
    print(f"  s histogram: {dict(s_pos)}")
    print(f"  r orbit-index histogram: {dict(r_pos)}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
