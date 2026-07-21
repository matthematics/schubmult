"""Verify the CLOSED-FORM exchange (correct convention) at n=5 AND n=6.

Moved case u NOT->_k w, u(k)!=w(k). A contributing up-neighbour is b>k with
w.inv - (w t_{kb}).inv in {1, 1-2(b-k)}  (i.e. w t_{kb} <| w) and u ->_{k-1} w t_{kb}.

Claims (expect 0 fails), whenever >=1 contributing up-neighbour exists:
 (E1) exactly one contributing up-neighbour, at b = B = max C ; and pi(k) == B.
 (E2) r := pi^{-1}(k) satisfies r < k, r in C, and t_{rk} is the UNIQUE contributing
      down-neighbour (a<k with w.inv-(w t_{ak}).inv in {1,1-2(k-a)} and u->_{k-1} w t_{ak}).
 (E3) P_{k-1}(u, w t_{kB}) == P_{k-1}(u, w t_{rk}) == P_{k-1}(u,w) ; equal n_{k-1}.
 (E4) kappa_{kB} q_{k-1}(u, w t_{kB}) == kappa_{rk} q_{k-1}(u, w t_{rk})  [expand].
 (E5) up split of C = (B fixed) + (cycle top k) ; down split = (k fixed) + (cycle top B),
      and the two nontrivial cycles have the SAME non-top element set {a_2,...,a_p}.
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
def Pset(u, w, m): return frozenset(val(u, i) for i in range(1, m+1) if val(u, i) == val(w, i))
def cyc_through(perm, x):
    c = [x]; cur = perm[x-1]
    while cur != x and cur not in c: c.append(cur); cur = perm[cur-1]
    return c
def cyc_containing(perm, x, N):
    c = [x]; cur = perm[x-1]
    while cur != x: c.append(cur); cur = perm[cur-1]
    return set(c)
def nontrivial_cycle_of(perm, x, N):
    if perm[x-1] == x: return None
    return frozenset(cyc_containing(perm, x, N))
def top_of_cycle(perm, x, N):
    return max(cyc_containing(perm, x, N))
def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    cases = 0; fE1 = fE2 = fE3 = fE4 = fE5 = 0
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
                        kap = S.One if (w.inv - vv.inv) == 1 else q_ab(k, b)
                        ups.append((b, kap, vv))
                if not ups: continue
                cases += 1
                C = cyc_through(pi, k); B = max(C)
                # E1
                if not (len(ups) == 1 and ups[0][0] == B and pi[k-1] == B):
                    fE1 += 1
                # E2
                downs = []
                for a in range(1, k):
                    vv = w.swap(a-1, k-1)
                    if (w.inv - vv.inv) in (1, 1-2*(k-a)) and vv in Rkm1:
                        kap = S.One if (w.inv - vv.inv) == 1 else q_ab(a, k)
                        downs.append((a, kap, vv))
                r = (~pi)[k-1]  # pi^{-1}(k)
                if not (len(downs) == 1 and downs[0][0] == r and r < k and r in C):
                    fE2 += 1
                if not downs:
                    continue
                b, kap_up, v_up = ups[0]; a, kap_dn, v_dn = downs[0]
                # E3
                if not (Pset(u, v_up, k-1) == Pset(u, v_dn, k-1) == Pset(u, w, k-1)):
                    fE3 += 1
                # E4
                if expand(kap_up*Rkm1[v_up]) != expand(kap_dn*Rkm1[v_dn]):
                    fE4 += 1
                # E5
                pu = pi.swap(k-1, B-1)   # pi t_{kB}
                pd = pi.swap(r-1, k-1)   # pi t_{rk}
                Dk = nontrivial_cycle_of(pu, k, N)   # cycle containing k after up
                DB = nontrivial_cycle_of(pd, B, N)   # cycle containing B after down
                okE5 = (pu[B-1] == B and pd[k-1] == k    # B fixed after up, k fixed after down
                        and Dk is not None and DB is not None
                        and top_of_cycle(pu, k, N) == k and top_of_cycle(pd, B, N) == B
                        and (Dk - {k}) == (DB - {B}))
                if not okE5:
                    fE5 += 1
    print(f"n={n}: cases with a contributing up-neighbour = {cases}")
    print(f"  (E1) unique up at b=B and pi(k)=B          : fails = {fE1}")
    print(f"  (E2) unique down at r=pi^{{-1}}(k), r<k, in C : fails = {fE2}")
    print(f"  (E3) P_{{k-1}} up==down==w                    : fails = {fE3}")
    print(f"  (E4) kappa*q weights match                 : fails = {fE4}")
    print(f"  (E5) split shapes (B fix|top k) / (k fix|top B), same lower set : fails = {fE5}")
if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
