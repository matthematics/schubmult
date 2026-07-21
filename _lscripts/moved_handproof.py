"""Verify the hand-proof sub-claims for the up-down exchange, and count down-neighbours.

Moved case u NOT->_k w, u(k)!=w(k), pi=u^{-1}w, C=cycle thru k, B=max C.
Consider (u,w,k) admitting a contributing UP-neighbour: some b>k with w t_{kb} <| w (lup) and
u ->_{k-1} w t_{kb}. Verify (expect 0 fails, n=5 and n=6):

 (H1) b == B == max C ; C has exactly two indices >= k, namely k and B ; pi(k)==B.
 (H2) u NOT ->_{k-1} w  (so DIAG,STRAIGHT vanish) ; w(k)==u(B).
 (H3) P_{k-1}(u, w t_{kB}) == P_{k-1}(u,w).
 (H4) for EVERY valid down-neighbour t_{ak} (a<k, w t_{ak} <| w, u->_{k-1}): a in C and
      P_{k-1}(u, w t_{ak}) == P_{k-1}(u,w)  and  n_{k-1} equal.
 (COUNT) number of valid down-neighbours, split by whether u->_{k-2}w (quantum present).
         Claim: #down == 1 + [u->_{k-2}w]. And exactly one down has q-weight == the up q-weight
         (kappa_kB q_{k-1}(u,v_up)); if u->_{k-2}w exactly one down has weight q_{k-1}q_{k-2}(u,w).
"""
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
                results.setdefault(nperm, nqw)
                stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results

def Pset(u, w, m): return frozenset(val(u, i) for i in range(1, m+1) if val(u, i) == val(w, i))
def npart(u, w, m): return m - len(Pset(u, w, m))
def cyc_through(perm, x, N):
    c = [x]; cur = perm[x-1]
    while cur != x and cur not in c: c.append(cur); cur = perm[cur-1]
    return c

def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    cases = 0
    fH1 = fH2 = fH3 = fH4 = 0
    down_count = Counter(); weightmatch_fail = 0
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rkm1 = enumerate_pieri(u, k-1, N); Rk = enumerate_pieri(u, k, N); Rkm2 = enumerate_pieri(u, k-2, N)
            cands = set(Rkm1) | set(Rk) | set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cands.add(v0.swap(k-1, pp-1))
            for w in cands:
                w = Permutation(w)
                if val(u, k) == val(w, k): continue
                if w in Rk: continue
                pi = (~u)*w
                # find contributing up-neighbours
                ups = []
                for b in range(k+1, N+1):
                    vv = w.swap(k-1, b-1); dl = w.inv - vv.inv
                    if dl == 1: kap = S.One
                    elif dl == 1-2*(b-k): kap = q_ab(k, b)
                    else: continue
                    if vv in Rkm1:
                        ups.append((b, kap, vv))
                if not ups:
                    continue
                cases += 1
                C = cyc_through(pi, k, N); B = max(C)
                geqk = [x for x in C if x >= k]
                # H1
                if not (len(ups) == 1 and ups[0][0] == B and set(geqk) == {k, B} and B > k and pi[k-1] == B):
                    fH1 += 1
                # H2
                if (w in Rkm1) or (val(w, k) != val(u, B)):
                    fH2 += 1
                b, kap_up, v_up = ups[0]
                # H3
                if Pset(u, v_up, k-1) != Pset(u, w, k-1):
                    fH3 += 1
                # down-neighbours
                downs = []
                for a in range(1, k):
                    vv = w.swap(a-1, k-1); dl = w.inv - vv.inv
                    if dl == 1: kap = S.One
                    elif dl == 1-2*(k-a): kap = q_ab(a, k)
                    else: continue
                    if vv in Rkm1:
                        downs.append((a, kap, vv))
                for a, kap, vv in downs:
                    if a not in C or Pset(u, vv, k-1) != Pset(u, w, k-1) or npart(u, vv, k-1) != npart(u, w, k-1):
                        fH4 += 1
                quantum = w in Rkm2
                down_count[(len(downs), quantum)] += 1
                # weight match: exactly one down has weight == up weight
                w_up = expand(kap_up*Rkm1[v_up])
                matches = [a for a, kap, vv in downs if expand(kap*Rkm1[vv]) == w_up]
                if len(matches) != 1:
                    weightmatch_fail += 1
    print(f"n={n} cases with an up-neighbour = {cases}")
    print(f"  (H1) b=B, C two indices {{k,B}}, pi(k)=B : fails = {fH1}")
    print(f"  (H2) u NOT->_{{k-1}}w & w(k)=u(B)         : fails = {fH2}")
    print(f"  (H3) P_{{k-1}}(u,wt_kB)=P_{{k-1}}(u,w)     : fails = {fH3}")
    print(f"  (H4) every down: a in C, P&n match       : fails = {fH4}")
    print(f"  (WT) not exactly one weight-matched down : fails = {weightmatch_fail}")
    print(f"  down-count histogram (ndown, u->_k-2 w)  : {dict(down_count)}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
