"""Does nabla (delete first cross in reading order) commute with D^p?  Compare nabla(D^p T) vs D^p(nabla T)."""
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

def nabla(X):
    rows = [list(r) for r in X]
    for r in rows:
        if r:
            r.pop(0)  # rows are stored decreasing; first in reading order = largest letter
            break
    return RCGraph([tuple(r) for r in rows])

def key(X):
    return tuple(tuple(r) for r in X)

stats = Counter()
bad = []
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(2, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                lhs = nabla(Dp(T, p))
                rhs = Dp(nabla(T), p)
                ok = key(lhs) == key(rhs)
                stats[ok] += 1
                if not ok and len(bad) < 4:
                    bad.append((rows, p, key(lhs), key(rhs)))
print(stats)
for b in bad:
    print(b)
