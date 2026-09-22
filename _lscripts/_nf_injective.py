"""Is (Q, seq, std(w(1..p))) injective on normal forms (height p RC graphs with maxd<=p)?"""
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

def std(vals):
    s = sorted(vals)
    return tuple(s.index(v) + 1 for v in vals)

tbl = defaultdict(set)
for N in range(2, 9):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        p = maxd(w)
        if p > 5:
            continue
        a = list(w) + list(range(len(w) + 1, N + 2))
        for A in RCGraph.all_rc_graphs(w, p):
            word, seq = A.as_reduced_compatible()
            key = (p, tuple(seq), eg_Q(tuple(reversed(word))), std(a[:p]))
            tbl[key].add(tuple(word))
amb = {k: v for k, v in tbl.items() if len(v) > 1}
print("keys", len(tbl), "ambiguous", len(amb))
for k, v in list(amb.items())[:6]:
    print(k, v)
