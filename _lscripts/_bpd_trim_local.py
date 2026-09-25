"""(a) Phi(shift-up T) == 1 (+) Phi(T)  (shift by (1,1) with straight pipe 1).
(b) For bounded X of height n: Phi(trim X)_{<= n-2} is determined by Phi(X)_{<= n-1}."""
import sys
from collections import Counter, defaultdict
from itertools import permutations

import numpy as np

from schubmult import BPD, Permutation, RCGraph

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
stats = Counter()


def phi(R, N):
    return BPD.from_rc_graph(R.resize(N)).resize(N)


for N0 in range(2, NMAX + 1):
    for a in permutations(range(1, N0 + 1)):
        w = Permutation(list(a))
        if w.inv == 0:
            continue
        d = len(w.trimcode)
        N = N0 + 2
        for n in range(max(d, 2), N0 + 1):
            groups = defaultdict(set)
            for X in RCGraph.all_rc_graphs(w, n):
                B = phi(X, N)
                # (a) shift
                up = RCGraph([()] + [tuple(l + 1 for l in r) for r in X])
                Bup = phi(up, N + 1)
                shifted = np.full((N + 1, N + 1), -1, dtype=object)
                shifted[1:, 1:] = B._grid
                ok = np.array_equal(Bup._grid[1:, 1:], B._grid)
                stats[("shift", ok)] += 1
                # (b)
                T = RCGraph([tuple(l - 1 for l in r) for r in list(X)[1:]])
                BT = phi(T, N)
                groups[B._grid[: n - 1].tobytes()].add(BT._grid[: n - 2].tobytes())
            for k, v in groups.items():
                stats[("trim local", len(v) == 1)] += 1
for k, v in sorted(stats.items(), key=str):
    print(v, k)
