"""For bounded RC graphs R of height n (maxd(perm) <= n): is phi^{-1}(R) (N x N) equal to the canonical completion
of its top n rows (resize(n).resize(N))?  Also: among BPDs of w with maxd(w) <= n and no blanks below row n, are the
top n rows injective?"""
import sys
from collections import Counter, defaultdict
from itertools import permutations

import numpy as np

from schubmult import BPD, Permutation, RCGraph

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
stats = Counter()
for N0 in range(2, NMAX + 1):
    for a in permutations(range(1, N0 + 1)):
        w = Permutation(list(a))
        if w.inv == 0:
            continue
        d = len(w.trimcode)
        N = N0 + 1
        for n in range(d, N0 + 1):
            tops = defaultdict(set)
            for R in RCGraph.all_rc_graphs(w, n):
                B = BPD.from_rc_graph(R.resize(N)).resize(N)
                assert B.perm == w
                C = B.resize(n).resize(N)
                stats[("canonical completion", np.array_equal(B._grid, C._grid))] += 1
                tops[B._grid[:n].tobytes()].add(B._grid.tobytes())
            for k, v in tops.items():
                stats[("top rows injective", len(v) == 1)] += 1
for k, v in sorted(stats.items(), key=str):
    print(v, k)
