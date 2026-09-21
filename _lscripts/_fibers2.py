"""Fibers of D^3 over R1=(2,3) and R2=(1,3) with seq (1,3): print T, w, low pattern, Inv_ll."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

def arr(w, n=9):
    a = list(w); a += list(range(len(a) + 1, n + 1)); return tuple(a[:n])

p = 3
fib = {(2, 3): [], (1, 3): []}
for N in range(2, 10):
    for a in permutations(range(1, N + 1)):
        w = Permutation(list(a))
        if w.inv != 2:
            continue
        md = maxd(w)
        for T in RCGraph.all_rc_graphs(w, max(md, p)):
            rows = tuple(tuple(r) for r in T)
            if any(len(rows[i]) for i in range(p, len(rows))):
                continue
            word, seq = T.as_reduced_compatible()
            if tuple(seq) != (1, 3):
                continue
            D = tuple(Dp(T, p).as_reduced_compatible()[0])
            if D in fib:
                aw = arr(w)
                fib[D].append((tuple(word), aw[:5]))
for k, v in fib.items():
    print("NF", k)
    for x in sorted(set(v)):
        print("   ", x)
