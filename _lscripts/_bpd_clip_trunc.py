"""Test BPD statements needed for Lemma clpfunc via Gao-Huang.

phi^{-1} = BPD.from_rc_graph, truncation = BPD.resize(p).
(Z)  zero_out_last_row(R) corresponds to BPD truncation by one row, for R with empty last row.
(T)  from_rc_graph(R).resize(p) == from_rc_graph(R_{<=p} padded).resize(p)   (top rows local)
(Z') from_rc_graph(clip^p(R)) == from_rc_graph(R).resize(p)
"""
import sys
from itertools import permutations
from collections import Counter
from schubmult import Permutation, RCGraph, BPD

def maxd(w):
    return len(w.trimcode)

def clip(R, p):
    rows = [tuple(r) for r in R][:p]
    T = RCGraph(rows)
    h = max(p, maxd(T.perm))
    T = T.resize(h)
    while len(T) > p:
        T = T.zero_out_last_row()
    return T

stats = Counter(); ex = {}
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
for N in range(2, NMAX + 1):
    for a in permutations(range(1, N + 1)):
        w = Permutation(list(a))
        if w.inv == 0:
            continue
        d = maxd(w)
        for n in range(d, N + 1):
            for R in RCGraph.all_rc_graphs(w, n):
                B = BPD.from_rc_graph(R)
                # (Z)
                if len(R[-1]) == 0 and n >= 2:
                    Z = R.zero_out_last_row()
                    ok = BPD.from_rc_graph(Z) == B.resize(n - 1)
                    stats[("Z", ok)] += 1
                    if not ok and "Z" not in ex:
                        ex["Z"] = ([tuple(r) for r in R], [tuple(r) for r in Z])
                for p in range(1, n):
                    rows = [tuple(r) for r in R][:p]
                    T = RCGraph(rows); h = max(p, maxd(T.perm)); T = T.resize(h)
                    okT = BPD.from_rc_graph(T).resize(p) == B.resize(p)
                    stats[("T", okT)] += 1
                    if not okT and "T" not in ex:
                        ex["T"] = ([tuple(r) for r in R], p)
                    C = clip(R, p)
                    okZp = BPD.from_rc_graph(C) == B.resize(p)
                    stats[("Z'", okZp)] += 1
                    if not okZp and "Z'" not in ex:
                        ex["Z'"] = ([tuple(r) for r in R], p, [tuple(r) for r in C])
for k, v in sorted(stats.items()):
    print(v, k)
for k, v in ex.items():
    print("counterexample", k, v)
