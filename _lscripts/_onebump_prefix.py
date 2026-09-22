"""Test (**): for Y with empty last row, maxd = height m, p < m:
either (lmap Y)_{<=p} == Y_{<=p}, or (lmap Y)_{<=p} == lmap(Y_{<=p}) (Y_{<=p} at height maxd(Y_{<=p}) > p)."""
from itertools import permutations
from schubmult import Permutation, RCGraph
from collections import Counter

def maxd(w):
    return len(w.trimcode)

def top(X, p):
    return RCGraph([tuple(r) for r in X][:p])

stats = Counter()
bad = []
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        m = maxd(w)
        if m < 3:
            continue
        for Y in RCGraph.all_rc_graphs(w, m):
            rows = [tuple(r) for r in Y]
            if len(rows[-1]) != 0:
                continue
            Y1 = Y.little_bump_desc()
            for p in range(2, m):
                T = top(Y, p)
                T1 = top(Y1, p)
                if T == T1:
                    stats["unchanged"] += 1
                    continue
                h = maxd(T.perm)
                if h <= p:
                    stats["FAIL: top bounded but changed"] += 1
                    if len(bad) < 5:
                        bad.append(("bounded", rows, p, [tuple(r) for r in Y1]))
                    continue
                Tb = T.resize(h).little_bump_desc()
                ok = top(Tb, p) == T1
                stats["bump-match" if ok else "FAIL: bump mismatch"] += 1
                if not ok and len(bad) < 5:
                    bad.append(("mismatch", rows, p, [tuple(r) for r in Y1], [tuple(r) for r in Tb]))
print(stats)
for b in bad:
    print(b)
