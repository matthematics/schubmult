"""Test the 'iterated exchange' reduction: for X = T u ^p S (height p+q), i a right descent of v = w_S,
remove the unique crossing of ^pS with root (p+i, p+i+1) (exchange property) to get X' = T u ^p S^{(i)}.
Is clip^p(X) == clip^p(X') (when X' is bounded at height p+q)?  Also test the recursive clip vs D^p on X'."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd
from collections import Counter

def one_row_clip(X):
    rows = [tuple(r) for r in X]; n = len(rows)
    top = RCGraph(rows[:n - 1]); h = max(n - 1, maxd(top.perm)); top = top.resize(h)
    while len(top) > n - 1:
        top = top.zero_out_last_row()
    return top

def clip(X, p):
    while len(X) > p:
        X = one_row_clip(X)
    return X

def Dp(X, p):
    rows = [tuple(r) for r in X][:p]
    T = RCGraph(rows); h = max(p, maxd(T.perm)); T = T.resize(h)
    while len(T) > p:
        T = T.zero_out_last_row()
    return T

stats = Counter(); ex = []
for N in range(2, 8):
    for a in permutations(range(1, N + 1)):
        w = Permutation(list(a))
        if w.inv > 6 or w.inv == 0:
            continue
        n = maxd(w)
        for X in RCGraph.all_rc_graphs(w, n):
            rows = [tuple(r) for r in X]
            for p in range(1, n):
                q = n - p
                Srows = [tuple(x - p for x in r) for r in rows[p:]]
                if all(len(r) == 0 for r in Srows):
                    continue
                S = RCGraph(Srows); v = S.perm
                for d in v.descents():  # 0-indexed descents
                    i = d + 1
                    S2, row = S.exchange_property(i, return_row=True)
                    # rebuild X' = T u ^p S2
                    rows2 = rows[:p] + [tuple(x + p for x in r) for r in S2]
                    X2 = RCGraph(rows2)
                    if maxd(X2.perm) > n:
                        stats["X' not bounded"] += 1
                        continue
                    if maxd(S2.perm) > q:
                        stats["S' not bounded (X' bounded)"] += 1
                    c1 = clip(X, p); c2 = clip(X2, p)
                    ok = c1 == c2
                    stats[("clip(X)==clip(X')", ok)] += 1
                    if not ok and len(ex) < 5:
                        ex.append((rows, p, i, rows2, [tuple(r) for r in c1], [tuple(r) for r in c2]))
print(stats)
for e in ex:
    print("  ", e)
