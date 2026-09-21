"""Tabulate w' = wof(D^p T) vs T for small T (rows<=p), to guess the low values w'(1..p)."""
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

out = []
for N in range(2, 7):
    for a in permutations(range(1, N + 1)):
        w = Permutation(list(a))
        if w.inv > 4 or w.inv == 0:
            continue
        md = maxd(w)
        for p in (2,):
            if md <= p:
                continue
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                D = Dp(T, p)
                out.append((rows[:p], arr(w), arr(D.perm), [tuple(r) for r in D]))
out.sort(key=lambda r: (len(r[0][0]) + len(r[0][1]), r[0]))
for r in out[:70]:
    print(f"T={r[0]} w={r[1]} -> D={r[3]} w'={r[2]}")
