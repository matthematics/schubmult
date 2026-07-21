"""Universal counting statement for pure-reflection groups (both B=k and B>k).

Moved case: u NOT->_k w, u(k)!=w(k). For each q-weight q^M, let G_M be the summands of the
recursion F of that weight. Classify: DIAG/STRAIGHT (weight q_{k-1}(u,w)); QUANT (weight
q_{k-1}q_{k-2}(u,w)); up-Monk t_{kb} b>k; down-Monk t_{ak} a<k.

CLAIM (verify 0-fail): every q-weight group G_M containing NEITHER diag/straight NOR quant
(a "pure reflection group") has #up_M == #down_M. Within a group all weights equal q^M by
construction, and (P_{k-1},n_{k-1}) common, so #up==#down  => group sums to 0.

Also report: do up-Monk neighbours EVER appear in a group that contains diag/straight or quant?
(Expect NO: those groups are handled by Case (i) / reflectioncancel with a unique down partner.)
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

def Pn(u, w, kk):
    P = frozenset(val(u, i) for i in range(1, kk+1) if val(u, i) == val(w, i))
    return P, kk-len(P)

def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    fail_count = 0; fail_pn = 0; up_in_special = 0; pure_groups = 0; moved = 0
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rkm1 = enumerate_pieri(u, k-1, N)
            Rkm2 = enumerate_pieri(u, k-2, N)
            Rk = enumerate_pieri(u, k, N)
            cset = set(Rkm1)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cset.add(v0.swap(k-1, pp-1))
            for w in cset:
                w = Permutation(w)
                if val(u, k) == val(w, k): continue
                if w in Rk: continue
                moved += 1
                # weights of diag/straight and quant
                specialw = set()
                if w in Rkm1: specialw.add(expand(Rkm1[w]))              # diag & straight
                if w in Rkm2: specialw.add(expand(q_gs[k-1]*q_gs[k-2]*Rkm2[w]))  # quant
                # gather Monk neighbours grouped by weight
                groups = {}   # weight -> {'up':[r], 'down':[r], 'pn':set()}
                for r in range(1, N+1):
                    if r == k: continue
                    vv = w.swap(k-1, r-1); dl = w.inv - vv.inv
                    lo, hi = min(r, k), max(r, k)
                    if dl == 1: kap = S.One
                    elif dl == 1-2*(hi-lo): kap = q_ab(lo, hi)
                    else: continue
                    if vv not in Rkm1: continue
                    wt = expand(kap*Rkm1[vv])
                    g = groups.setdefault(wt, {'up': [], 'down': [], 'pn': set()})
                    (g['up'] if r > k else g['down']).append(r)
                    g['pn'].add(Pn(u, vv, k-1))
                for wt, g in groups.items():
                    if wt in specialw:
                        if g['up']:
                            up_in_special += 1     # an up-neighbour in a diag/quant group
                        continue
                    # pure reflection group
                    pure_groups += 1
                    if len(g['up']) != len(g['down']):
                        fail_count += 1
                    if len(g['pn']) != 1:
                        fail_pn += 1
    print(f"n={n} moved cases={moved} pure-reflection groups={pure_groups}")
    print(f"  #up != #down in a pure group : fails = {fail_count}")
    print(f"  (P_,n_) not common in a pure group : fails = {fail_pn}")
    print(f"  up-neighbour inside a diag/straight/quant group : {up_in_special}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
