"""Test C2: R == zeromap^{N-n}( R_{<=p} star R_{>p} ), where star shifts rows p+1..n far right."""
from itertools import permutations
from schubmult import Permutation, RCGraph
from collections import Counter

def maxd(w):
    return len(w.trimcode)

def zero_to(X, p):
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

def star(rows, p, M):
    # rows: list of tuples of letters; shift letters in rows > p by M
    return [tuple(r) if i < p else tuple(x + M for x in r) for i, r in enumerate(rows)]

stats = Counter()
bad = []
for N in range(2, 7):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 6 or w.inv == 0:
            continue
        md = maxd(w)
        for n in range(md, md + 2):
            for R in RCGraph.all_rc_graphs(w, n):
                rows = [tuple(r) for r in R]
                maxletter = max([x for r in rows for x in r])
                for p in [n - 1]:
                    if all(len(rows[i]) == 0 for i in range(p, n)):
                        continue
                    M = maxletter + 2
                    srows = star(rows, p, M)
                    X = RCGraph(srows)
                    h = max(n, maxd(X.perm))
                    X = X.resize(h)
                    Z = zero_to(X, n)
                    ok = [tuple(r) for r in Z] == rows
                    stats[ok] += 1
                    if not ok and len(bad) < 5:
                        bad.append((rows, p, [tuple(r) for r in Z]))
print(stats)
for b in bad:
    print(b)
