"""Single-step forward locality of Gao-Huang's pop (nabla):
for a BPD B whose top blank row is r <= p, compare pop(B) with pop(C), C = completion of B_{<=p}.
Check (a) popped letter equal, (b) top p rows of nabla(B) and nabla(C) equal."""
import sys
from collections import Counter

import numpy as np

from schubmult import Permutation
from schubmult.combinatorics.bpd import BPD, TileType

n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
stats = Counter()
bad = []
for w in Permutation.all_permutations(n):
    if w.inv == 0:
        continue
    for B in BPD.all_bpds(w):
        B = B.resize(n)
        blanks = np.argwhere(B._grid == TileType.BLANK)
        r = int(blanks[0, 0]) + 1
        DB, (aB, rB) = B.pop_op()
        for p in range(r, n):
            C = B.resize(p).resize(n)
            assert np.array_equal(C._grid[:p], B._grid[:p])
            DC, (aC, rC) = C.pop_op()
            ok_a = (aB, rB) == (aC, rC)
            ok_top = np.array_equal(DB.resize(n)._grid[:p], DC.resize(n)._grid[:p])
            stats[("letter", ok_a)] += 1
            stats[("top rows", ok_top)] += 1
            if not (ok_a and ok_top) and len(bad) < 3:
                bad.append((w, p, B, C, (aB, rB), (aC, rC)))
for k, v in sorted(stats.items(), key=str):
    print(v, k)
for w, p, B, C, pb, pc in bad:
    print("w", w, "p", p, "pop(B)", pb, "pop(C)", pc)
    print(B)
    print(C)
