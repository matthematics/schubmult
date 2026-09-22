"""Can rebounding always be avoided? For Y of height n (n >= maxd(w_Y)) and p<n with rows>p nonempty:
does there exist a crossing in rows > p at a simple root (j,j+1) such that maxd(w s_j) <= n?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd, rroots
from collections import Counter

stats = Counter(); ex = []
for N in range(2, 8):
    for a in permutations(range(1, N + 1)):
        w = Permutation(list(a))
        if w.inv > 6 or w.inv == 0:
            continue
        m = maxd(w)
        for n in (m, m + 1):
            for Y in RCGraph.all_rc_graphs(w, n):
                rows = [tuple(r) for r in Y]
                word, seq = Y.as_reduced_compatible(); R = rroots(word)
                for p in range(1, n):
                    if all(len(rows[i]) == 0 for i in range(p, n)):
                        continue
                    cands = []
                    for k in range(len(word)):
                        if seq[k] <= p:
                            continue
                        a_, b_ = sorted(R[k])
                        if b_ == a_ + 1:
                            w2 = w * Permutation.ref_product(a_)
                            cands.append((a_, maxd(w2) <= n))
                    ok = any(c[1] for c in cands)
                    stats[("some non-rebounding removal exists", ok)] += 1
                    if not ok and len(ex) < 6:
                        ex.append((rows, p, n, list(w), cands))
print(stats)
for e in ex:
    print("  ", e)
