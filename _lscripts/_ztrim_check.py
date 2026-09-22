"""Sanity check of lemma ztrim: trim(zeromap Y) == zeromap(trim Y) for Y with empty last row, height m>=2."""
from itertools import permutations
from schubmult import Permutation, RCGraph
from collections import Counter

def maxd(w):
    return len(w.trimcode)

def trim(X):
    rows = [tuple(a - 1 for a in r) for r in list(X)[1:]]
    return RCGraph(rows)

def zero(X):
    # X has empty last row; if maxd < height, plain deletion
    return X.zero_out_last_row()

stats = Counter()
bad = []
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for m in (md, md + 1):
            if m < 2:
                continue
            for Y in RCGraph.all_rc_graphs(w, m):
                rows = [tuple(r) for r in Y]
                if len(rows[-1]) != 0:
                    continue
                lhs = trim(zero(Y))
                T = trim(Y)
                rhs = zero(T.resize(m - 1))
                ok = [tuple(r) for r in lhs] == [tuple(r) for r in rhs]
                stats[ok] += 1
                if not ok and len(bad) < 5:
                    bad.append((rows, [tuple(r) for r in lhs], [tuple(r) for r in rhs]))
print(stats)
for b in bad:
    print(b)
