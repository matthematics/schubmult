"""Test: D^p(little_bump_desc(Y)) == D^p(Y) for Y with empty last row and maxd = height."""
from itertools import permutations
from schubmult import Permutation, RCGraph
from collections import Counter

def maxd(w):
    return len(w.trimcode)

def Dp(X, p):
    rows = [tuple(r) for r in X][:p]
    top = RCGraph(rows)
    h = max(p, maxd(top.perm))
    top = top.resize(h)
    while len(top) > p:
        top = top.zero_out_last_row()
    return top

stats = Counter()
bad = []
for N in range(2, 7):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 6 or w.inv == 0:
            continue
        m = maxd(w)
        if m < 3:
            continue
        for Y in RCGraph.all_rc_graphs(w, m):  # height m, need empty last row
            rows = [tuple(r) for r in Y]
            if len(rows[-1]) != 0:
                continue
            Y1 = Y.little_bump_desc()
            r1 = [tuple(r) for r in Y1]
            assert len(r1) == m and len(r1[-1]) == 0, (rows, r1)
            for p in range(2, m):
                ok = Dp(Y1, p) == Dp(Y, p)
                stats[ok] += 1
                if not ok and len(bad) < 5:
                    bad.append((rows, p, r1, [tuple(r) for r in Dp(Y, p)], [tuple(r) for r in Dp(Y1, p)]))
print(stats)
for b in bad:
    print(b)
