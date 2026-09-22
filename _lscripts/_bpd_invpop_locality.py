"""Backward locality of the inverse pop (BL1):
if E, E' agree in rows <= q and D = nabla^{-1}_{(a,r)}(E), D' = nabla^{-1}_{(a,r)}(E') with r <= q, then D, D' agree in rows <= q.
Tested by grouping all D by (top q rows of nabla D, pop(D))."""
import sys
from collections import Counter, defaultdict

import numpy as np

from schubmult import Permutation
from schubmult.combinatorics.bpd import BPD, TileType

n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
stats = Counter()
bad = []
for q in range(1, n):
    groups = defaultdict(set)
    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        for D in BPD.all_bpds(w):
            D = D.resize(n)
            E, (a, r) = D.pop_op()
            if r > q:
                continue
            E = E.resize(n)
            key = (E._grid[:q].tobytes(), a, r)
            groups[key].add(D._grid[:q].tobytes())
    for key, s in groups.items():
        stats[(q, "BL1", len(s) == 1)] += 1
        if len(s) > 1 and len(bad) < 3:
            bad.append((q, key, s))
for k, v in sorted(stats.items(), key=str):
    print(v, k)
print("bad examples:", len(bad))
for q, key, s in bad:
    _, a, r = key
    print("q", q, "pop", (a, r))
    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        for D in BPD.all_bpds(w):
            D = D.resize(n)
            E, pr = D.pop_op()
            if pr == (a, r) and E.resize(n)._grid[:q].tobytes() == key[0]:
                print("D (perm", D.perm, ") -> E (perm", E.perm, ")")
                print(D)
                print(E.resize(n))
