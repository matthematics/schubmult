"""Tests:
(M) Are chains of general p-high bumps on RC graphs (rows<=p) strictly backward-monotone in position?
(U) Uniqueness: A, B in RC_p (height p, maxd<=p), same weight, same rows 2..p, D^{p-1}(A)==D^{p-1}(B)  =>  A==B ?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, maxd
from collections import Counter, defaultdict

def Dp(X, p):
    rows = [tuple(r) for r in X][:p]
    top = RCGraph(rows)
    h = max(p, maxd(top.perm))
    top = top.resize(h)
    while len(top) > p:
        top = top.zero_out_last_row()
    return top

# (M)
statsM = Counter()
badM = []
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(1, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible()
                for q in range(len(word)):
                    a, b = rroot(word, q)
                    if a <= p:
                        continue
                    nw, chain = ltb_down(word, q, return_chain=True)
                    if nw is None:
                        statsM["degenerate"] += 1
                        continue
                    dec = all(chain[k] > chain[k + 1] for k in range(len(chain) - 1))
                    statsM[("monotone", dec)] += 1
                    if not dec and len(badM) < 5:
                        badM.append((rows, p, q, (a, b), chain))
print("(M)", statsM)
for b in badM:
    print(b)

# (U)
tbl = defaultdict(set)
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7:
            continue
        p = maxd(w)
        if p < 2:
            continue
        for A in RCGraph.all_rc_graphs(w, p):
            rows = [tuple(r) for r in A]
            key = (p, tuple(len(r) for r in rows), tuple(rows[1:]), tuple(tuple(r) for r in Dp(A, p - 1)))
            tbl[key].add(rows[0])
amb = {k: v for k, v in tbl.items() if len(v) > 1}
print("(U) keys", len(tbl), "ambiguous", len(amb))
for k, v in list(amb.items())[:5]:
    print(k, v)
