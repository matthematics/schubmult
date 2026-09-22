"""Compare little_bump_desc (paper's lmap) with the word-level Little bump at the position with right-root (m,m+1)."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, maxd
from collections import Counter

stats = Counter(); ex = []
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        m = maxd(w)
        if m < 2:
            continue
        for Y in RCGraph.all_rc_graphs(w, m):
            rows = [tuple(r) for r in Y]
            if len(rows[-1]) != 0:
                continue
            word, seq = Y.as_reduced_compatible()
            word = tuple(word)
            Y1 = Y.little_bump_desc()
            w1 = tuple(Y1.as_reduced_compatible()[0])
            qs = [q for q in range(len(word)) if rroot(word, q) == (m, m + 1)]
            assert len(qs) == 1
            q = qs[0]
            nw, ch = ltb_down(word, q, return_chain=True)
            if nw is None:
                stats["degenerate"] += 1; continue
            ok = nw == w1
            inc = any(ch[k] < ch[k + 1] for k in range(len(ch) - 1))
            stats[("lmap==ltb", ok, "chain has increase", inc)] += 1
            if not ok and len(ex) < 5:
                ex.append((rows, word, q, ch, nw, w1, [tuple(r) for r in Y1]))
print(stats)
for e in ex:
    print(e)
