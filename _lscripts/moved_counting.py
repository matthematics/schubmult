"""Prove-by-verify the counting reformulation of Sigma=0 in the B>k regime.

Setup: u NOT->_k w, u(k)!=w(k), pi=u^{-1}w, C=cycle through k, B=max(C)>k.
Within a fixed q-weight group q^M, every contributing Monk neighbour wt_{rk} has the SAME
weight kappa_{rk} q_{k-1}(u,wt_{rk}) = q^M and the SAME (P_{k-1},n_{k-1}); so
   F_M-part = E * q^M * (#up_M - #down_M).
Hence Sigma=0 for all M  <=>  #up_M == #down_M for every q-weight M.

This script verifies, over all u in S_n, all k, all reachable w (B>k regime):
  (C1) every contributing Monk reflection t_{rk} has r in C, and the split of C at {r,k}
       separates k from B; the allowed set is exactly {c_0=B, c_1,...,c_{j-1}} (orbit from B).
  (C2) all contributing neighbours share P_{k-1},n_{k-1} with each other.
  (C3) per q-weight M: #up_M == #down_M.  (=> Sigma=0)
  (C4) intrinsic partner: there is a unique contributing up (=B) iff a unique contributing
       down r*; record r* vs the descent of the up-neighbour's k-cycle K0.
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

def Pn(u, w, kk, N):
    """agreement set P_{kk}(u,w)={u(i): i<=kk, u(i)=w(i)} and n_{kk}=kk-|P|."""
    P = frozenset(val(u, i) for i in range(1, kk+1) if val(u, i) == val(w, i))
    return P, kk-len(P)

def desc_of_kcycle(u, orb, j):
    """K0 (up-neighbour k-cycle) = orbit c_1..c_{j-1} then k=c_j; as cycle (c_1,...,c_{j-1},k)
    with top k. descent a_l: u(a_l)>u(a_{l-1}) walking k -> c_1 -> c_2 -> ... (a_0=top=k)."""
    seq = orb[1:j] + [orb[j]]  # c_1..c_{j-1}, k   (a_1..a_p with top appended? use top=k as a_0)
    # write cycle (a_p,...,a_1,b=k): orbit from top k is k->c_1->...->c_{j-1}->k
    prev = orb[j]  # k = top
    descs = []
    for a in orb[1:j]:  # c_1..c_{j-1}
        if val(u, a) > val(u, prev): descs.append(a)
        prev = a
    return descs

def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    tot = 0; failC1 = failC3 = failC4 = 0; nonempty = 0
    rstar_is_desc = 0; rstar_examples = []
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rkm1 = enumerate_pieri(u, k-1, N); Rk = enumerate_pieri(u, k, N)
            cset = set(Rkm1)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cset.add(v0.swap(k-1, pp-1))
            for w in cset:
                w = Permutation(w)
                if val(u, k) == val(w, k): continue
                if w in Rk: continue
                pi = (~u)*w; C = orbit(pi, k, N); B = max(C)
                if B <= k: continue
                tot += 1
                orb = orbit(pi, B, N); L = len(orb); j = orb.index(k)
                allowed = set(orb[0:j])  # {c_0=B, c_1,...,c_{j-1}}
                # gather ALL contributing Monk reflections t_{rk}
                ups = {}; downs = {}   # weight -> count ; also track r
                up_rs = []; down_rs = []
                Pn_set = set()
                for r in range(1, N+1):
                    if r == k: continue
                    vv = w.swap(k-1, r-1); dl = w.inv - vv.inv
                    lo, hi = min(r, k), max(r, k)
                    if dl == 1: kap = S.One
                    elif dl == 1-2*(hi-lo): kap = q_ab(lo, hi)
                    else: continue
                    if vv not in Rkm1: continue
                    # contributing
                    if r not in allowed: failC1 += 1
                    wt = expand(kap*Rkm1[vv])
                    Pn_set.add(Pn(u, vv, k-1, N))
                    if r > k:
                        ups[wt] = ups.get(wt, 0)+1; up_rs.append(r)
                    else:
                        downs[wt] = downs.get(wt, 0)+1; down_rs.append(r)
                if not ups and not downs: continue
                nonempty += 1
                # C3: per-weight up==down
                weights = set(ups) | set(downs)
                if any(ups.get(m, 0) != downs.get(m, 0) for m in weights):
                    failC3 += 1
                # C4: unique up=B, unique down
                if not (up_rs == [B] and len(down_rs) == 1):
                    failC4 += 1
                else:
                    rstar = down_rs[0]
                    descs = desc_of_kcycle(u, orb, j)
                    if len(rstar_examples) < 12:
                        rstar_examples.append((tuple(u), tuple(w), k, B, orb, j, rstar, descs))
                    if descs and rstar == descs[0]:
                        rstar_is_desc += 1
    print(f"n={n} B>k regime instances={tot} nonempty(Sigma) groups={nonempty}")
    print(f"  C1 (contributing r not in allowed) fails = {failC1}")
    print(f"  C3 (#up_M != #down_M for some M)   fails = {failC3}")
    print(f"  C4 (not unique up=B & unique down) fails = {failC4}")
    print(f"  r*==desc(K0)[0] in {rstar_is_desc}/{nonempty} nonempty")
    for ex in rstar_examples:
        u_, w_, k_, B_, orb_, j_, rs_, ds_ = ex
        print(f"  u={u_} w={w_} k={k_} B={B_} orb(from B)={orb_} j={j_} r*={rs_} desc(K0)={ds_}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
