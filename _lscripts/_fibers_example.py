"""Print fibers of D^p for the two ambiguous normal forms (2,1,3) and (4,1,3), p=3, seq (1,1,3)."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd

def Dp(X, p):
    h = max(p, maxd(X.perm))
    X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

def arr(w, n=8):
    a = list(w); a += list(range(len(a) + 1, n + 1)); return tuple(a[:n])

p = 3
targets = {((2, 1), (), (3,)): [], ((4, 1), (), (3,)): []}
for N in range(2, 9):
    for a in permutations(range(1, N + 1)):
        w = Permutation(list(a))
        if w.inv != 3:
            continue
        md = maxd(w)
        for T in RCGraph.all_rc_graphs(w, max(md, p)):
            rows = tuple(tuple(r) for r in T)
            if any(len(rows[i]) for i in range(p, len(rows))):
                continue
            if tuple(len(r) for r in rows[:p]) != (2, 0, 1):
                continue
            D = tuple(tuple(r) for r in Dp(T, p))
            if D in targets:
                targets[D].append((rows[:p], arr(w)))
for k, v in targets.items():
    print("NF", k)
    for x in sorted(v):
        print("   ", x)
