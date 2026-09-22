"""Relation R = {(Y_{<=p}, (ZY)_{<=p})}. Test: for (T,T') in R and position k (in rows<=p) with both roots p-high covers,
let X = bump of Y at position k (if the root of k in Y is a cover of w_Y and chain stays in rows<=p).
Is (Z X)_{<=p} == beta'(T')  (i.e. R preserved by simultaneous bumps)?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, valid_rc, maxd
from collections import Counter

def perm_of(word):
    return Permutation.ref_product(*word) if word else Permutation([])
def arr(w, n=14):
    a = list(w); a += list(range(len(a) + 1, n + 1)); return a[:n]
def is_cover(w, a, b):
    ar = arr(w)
    if ar[a - 1] < ar[b - 1]: return False
    return not any(ar[b - 1] < ar[c - 1] < ar[a - 1] for c in range(a + 1, b))

stats = Counter(); ex = []
for N in range(2, 9):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0: continue
        m = maxd(w)
        for p in range(2, m - 1):
            for Y in RCGraph.all_rc_graphs(w, m):
                rows = [tuple(r) for r in Y]
                if len(rows[p]) == 0 or any(len(rows[i]) for i in range(p + 1, m)): continue
                ZY = Y.zero_out_last_row()
                wy, sy = Y.as_reduced_compatible(); wy = tuple(wy); sy = tuple(sy)
                T = RCGraph(rows[:p]); T2 = RCGraph([tuple(r) for r in ZY][:p])
                wa, sa = T.as_reduced_compatible(); wb, _ = T2.as_reduced_compatible()
                wa, wb, sa = tuple(wa), tuple(wb), tuple(sa)
                u, v = perm_of(wa), perm_of(wb)
                for k in range(len(wa)):
                    ra, rb = rroot(wa, k), rroot(wb, k)
                    if not (ra[0] > p and is_cover(u, *ra) and rb[0] > p and is_cover(v, *rb)): continue
                    nb, _ = ltb_down(wb, k, return_chain=True)
                    if nb is None or not valid_rc(nb, sa): continue
                    ry = rroot(wy, k)
                    covY = is_cover(w, *ry)
                    stats[("root of k in Y is cover of w_Y", covY)] += 1
                    if not covY: continue
                    ny, Cy = ltb_down(wy, k, return_chain=True)
                    if ny is None or not valid_rc(ny, sy): stats["degenerate on Y"] += 1; continue
                    stays = all(sy[c] <= p for c in Cy)
                    stats[("chain of Y-bump stays in rows<=p", stays)] += 1
                    X = RCGraph.from_reduced_compatible(list(ny), list(sy))
                    if len(X) < m: X = X.resize(m)
                    if maxd(X.perm) != m or len([tuple(r) for r in X][-1]) != 0:
                        stats[("X has maxd==m and empty last row", False)] += 1; continue
                    ZX = X.zero_out_last_row()
                    lhs = tuple(x for r in [tuple(r) for r in ZX][:p] for x in r)
                    ok = lhs == nb
                    stats[("(ZX)_{<=p} == beta'(T')", ok, "chain stays", stays)] += 1
                    if not ok and len(ex) < 4:
                        ex.append((rows[:p + 1], p, k, wa, wb, ny[:len(wa)], lhs, nb))
for k, v in sorted(stats.items(), key=lambda x: str(x[0])):
    print(v, k)
for e in ex: print("  ", e)
