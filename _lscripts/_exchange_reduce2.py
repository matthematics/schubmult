"""(E) with rebounding: Y in RC_n, x a crossing in rows > p whose root is simple (j,j+1) (so Y\\x reduced).
Y' = Y \\ x regarded at height max(n, maxd(w s_j)). Is clip^p(Y') == clip^p(Y)? Split by rebounded or not, and by row of x."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd, rroots
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

stats = Counter(); ex = []
for N in range(2, 8):
    for a in permutations(range(1, N + 1)):
        w = Permutation(list(a))
        if w.inv > 6 or w.inv == 0:
            continue
        n = maxd(w)
        for Y in RCGraph.all_rc_graphs(w, n):
            rows = [tuple(r) for r in Y]
            word, seq = Y.as_reduced_compatible()
            R = rroots(word)
            for p in range(1, n):
                cY = None
                for k in range(len(word)):
                    if seq[k] <= p:
                        continue
                    a_, b_ = sorted(R[k])
                    if b_ != a_ + 1:
                        continue
                    # remove crossing k
                    nw = list(word[:k]) + list(word[k + 1:]); ns = list(seq[:k]) + list(seq[k + 1:])
                    Yp = RCGraph.from_reduced_compatible(nw, ns).resize(n)
                    h2 = max(n, maxd(Yp.perm))
                    Yp = Yp.resize(h2)
                    if cY is None:
                        cY = clip(Y, p)
                    c2 = clip(Yp, p)
                    ok = (c2 == cY)
                    last_row = (seq[k] == n)
                    stats[("rebounded", h2 > n, "x in last row", last_row, "equal", ok)] += 1
                    if not ok and len(ex) < 5:
                        ex.append((rows, p, k, (a_, b_), [tuple(r) for r in Yp], [tuple(r) for r in cY], [tuple(r) for r in c2]))
for k, v in sorted(stats.items()):
    print(v, k)
for e in ex:
    print("  ", e)
