"""Which single-letter deletions commute with D^p?  For T (rows<=p), position k with T\\k reduced:
D^p(T\\k) == D^p(T)\\k ?  Classify by whether k is first-in-row / last-in-row / interior."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd
from schubmult.utils.perm_utils import is_reduced
from collections import Counter

def Dp_ws(word, seq, p):
    T = RCGraph.from_reduced_compatible(list(word), list(seq))
    h = max(p, maxd(T.perm))
    T = T.resize(h)
    while len(T) > p:
        T = T.zero_out_last_row()
    w2, s2 = T.as_reduced_compatible()
    return tuple(w2), tuple(s2)

stats = Counter()
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv < 2:
            continue
        md = maxd(w)
        for p in range(1, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible()
                word = tuple(word); seq = tuple(seq)
                D = Dp_ws(word, seq, p)
                for k in range(len(word)):
                    nw = word[:k] + word[k + 1:]; ns = seq[:k] + seq[k + 1:]
                    if not is_reduced(list(nw)):
                        continue
                    Dk = Dp_ws(nw, ns, p)
                    ok = Dk == (D[0][:k] + D[0][k + 1:], D[1][:k] + D[1][k + 1:])
                    first_in_row = (k == 0 or seq[k - 1] != seq[k])
                    last_in_row = (k == len(word) - 1 or seq[k + 1] != seq[k])
                    kind = ("k=0" if k == 0 else "first-in-row" if first_in_row else "last-in-row" if last_in_row else "interior")
                    stats[(kind, ok)] += 1
print(sorted(stats.items()))
