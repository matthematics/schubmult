"""Test: D^p(x.T) determined by (x, D^p(T))?  Group RC graphs T (rows<=p) by (x, D^p(nabla T)) and check D^p(T) unique."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd
from collections import defaultdict

def Dp_ws(word, seq, p):
    if not word:
        return ((), ())
    T = RCGraph.from_reduced_compatible(list(word), list(seq))
    h = max(p, maxd(T.perm))
    T = T.resize(h)
    while len(T) > p:
        T = T.zero_out_last_row()
    w2, s2 = T.as_reduced_compatible()
    return (tuple(w2), tuple(s2))

tbl = defaultdict(set)
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(1, md + 1):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible()
                word = tuple(word); seq = tuple(seq)
                D = Dp_ws(word, seq, p)
                D0 = Dp_ws(word[1:], seq[1:], p)
                tbl[(p, word[0], seq[0], D0)].add(D)
amb = {k: v for k, v in tbl.items() if len(v) > 1}
print("keys", len(tbl), "ambiguous", len(amb))
for k, v in list(amb.items())[:6]:
    print(k, v)
