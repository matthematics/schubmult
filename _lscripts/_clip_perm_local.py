from itertools import permutations
from schubmult import Permutation, uncode, RCGraph
from collections import defaultdict

def maxd(w):
    c = list(w.trimcode)
    return len(c)

def D(X, p):
    # X: RCGraph of height h with rows > p empty; zero down to height p
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

tbl = defaultdict(set)
count = 0
for N in range(2, 7):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 6:
            continue
        md = maxd(w)
        if md == 0:
            continue
        for h in range(md, md + 2):
            for X in RCGraph.all_rc_graphs(w, h):
                rows = [tuple(r) for r in X]
                for p in range(2, h):
                    if all(len(rows[i]) == 0 for i in range(p, h)):
                        Y = D(X, p)
                        key = (p, tuple(w), tuple(len(r) for r in rows[:p]))
                        tbl[key].add(tuple(Y.perm))
                        count += 1
bad = {k: v for k, v in tbl.items() if len(v) > 1}
print("checked", count, "keys", len(tbl), "ambiguous", len(bad))
for k, v in list(bad.items())[:5]:
    print(k, v)
