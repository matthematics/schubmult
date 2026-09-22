"""Find a small example where the first one-row clip of T u ^p S changes the remaining rows of S."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd

def one_row_clip(X):
    rows = [tuple(r) for r in X]
    n = len(rows)
    top = RCGraph(rows[:n - 1])
    h = max(n - 1, maxd(top.perm))
    top = top.resize(h)
    while len(top) > n - 1:
        top = top.zero_out_last_row()
    return top

p, q = 2, 2
found = 0
for N in range(3, 7):
    for a in permutations(range(1, N + 1)):
        w = Permutation(list(a))
        if w.inv > 5 or maxd(w) > p + q:
            continue
        for X in RCGraph.all_rc_graphs(w, p + q):
            rows = [tuple(r) for r in X]
            if len(rows[p]) == 0 or len(rows[p + 1]) == 0:
                continue
            X1 = one_row_clip(X)
            r1 = [tuple(r) for r in X1]
            if r1[p] != rows[p] or r1[:p] != rows[:p]:
                T = rows[:p]; S = [tuple(x - p for x in r) for r in rows[p:]]
                print(f"w={''.join(map(str,list(w)))}  T={T}  S={S}  X=T u ^2S={rows}   one-row clip -> {r1}   (S's remaining row was {rows[p]}, now {r1[p]}; T was {T}, now {r1[:p]})")
                found += 1
                if found >= 6:
                    raise SystemExit
