"""(T1) For an RC graph R (height n >= len(perm) so phi^{-1} is reliable) with x = a cross in the last nonempty row alpha,
compare phi^{-1}(R) and phi^{-1}(R \\ x) in rows < alpha.
Variants: x = leftmost cross of row alpha (T1), x = rightmost (T1r), x = any (T1any).
Also record whether R \\ x's BPD move is consistent. Heights are padded to N = max(n, len(perm)+1)."""
import sys
from collections import Counter
from itertools import permutations

import numpy as np

from schubmult import BPD, Permutation, RCGraph

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
stats = Counter()
ex = {}


def bpd_of(rows, N):
    R = RCGraph([tuple(r) for r in rows] + [()] * (N - len(rows)))
    return BPD.from_rc_graph(R).resize(N)


for N0 in range(2, NMAX + 1):
    for a in permutations(range(1, N0 + 1)):
        w = Permutation(list(a))
        if w.inv == 0:
            continue
        for R in RCGraph.all_rc_graphs(w, N0):
            rows = [tuple(r) for r in R]
            nonempty = [i for i, r in enumerate(rows) if len(r) > 0]
            alpha = nonempty[-1]  # 0-based index of last nonempty row
            if alpha == 0:
                continue
            N = N0 + 1
            B = bpd_of(rows, N)
            row = rows[alpha]  # letters, decreasing; leftmost cross = smallest letter = last entry
            for tag, letter in [("T1 leftmost", row[-1]), ("T1r rightmost", row[0])] + [("T1any", l) for l in row]:
                new_rows = list(rows)
                new_rows[alpha] = tuple(l for l in row if l != letter)
                B2 = bpd_of(new_rows, N)
                ok = np.array_equal(B._grid[:alpha], B2._grid[:alpha])
                stats[(tag, ok)] += 1
                if not ok and tag not in ex:
                    ex[tag] = (rows, alpha + 1, letter)
for k, v in sorted(stats.items()):
    print(v, k)
for k, v in ex.items():
    print("counterexample", k, v)
