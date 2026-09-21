"""(U') A,B height p RC graphs with maxd<=p, same weight, same Q (EG recording), same word after deleting first letter
=> A==B?  Also check the version without Q."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd
from collections import defaultdict

def eg_Q(word):
    P, Q = [], []
    for t, x in enumerate(word, 1):
        r = 0
        while True:
            if r == len(P):
                P.append([x]); Q.append([t]); break
            row = P[r]
            bigger = [k for k, y in enumerate(row) if y > x]
            if not bigger:
                row.append(x); Q[r].append(t); break
            k = bigger[0]; y = row[k]
            if y == x + 1 and k > 0 and row[k - 1] == x:
                x = y
            else:
                row[k] = x; x = y
            r += 1
    return tuple(tuple(r) for r in Q)

tbl = defaultdict(set); tbl0 = defaultdict(set)
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        p = maxd(w)
        for A in RCGraph.all_rc_graphs(w, p):
            word, seq = A.as_reduced_compatible()
            key0 = (p, tuple(seq), tuple(word[1:]))
            tbl0[key0].add(word[0])
            # EG insertion in HY convention inserts right-to-left; use reversed word for Q
            key = (p, tuple(seq), tuple(word[1:]), eg_Q(tuple(reversed(word))))
            tbl[key].add(word[0])
amb0 = {k: v for k, v in tbl0.items() if len(v) > 1}
amb = {k: v for k, v in tbl.items() if len(v) > 1}
print("without Q: keys", len(tbl0), "ambiguous", len(amb0))
print("with Q: keys", len(tbl), "ambiguous", len(amb))
for k, v in list(amb.items())[:5]:
    print(k, v)
