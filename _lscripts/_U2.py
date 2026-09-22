"""(U''): highest-weight A, B in RC_p (height p, maxd<=p) with trim A == trim B, same weight, same extwt  =>  A == B ?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations, combinations
from schubmult import Permutation, RCGraph
from schubmult.utils.perm_utils import is_reduced
from _lb_common import maxd
from collections import defaultdict, Counter

def is_hw(T, p):
    return all(T.raising_operator(i) is None for i in range(1, p))

def extwt(T, p):
    seen = {T}; stack = [T]; best = tuple(len(r) for r in T)[:p]
    while stack:
        X = stack.pop()
        wt = tuple(len(r) for r in X)[:p]
        if wt < best:
            best = wt
        for i in range(1, p):
            Y = X.lowering_operator(i)
            if Y is not None and Y not in seen:
                seen.add(Y); stack.append(Y)
    return best

def rc_graphs(p, L):
    """All RC graphs with p rows, L crossings, letters <= L+p, valid, with maxd(perm) <= p."""
    def rows_gen(i, remaining):
        if i > p:
            if remaining == 0:
                yield []
            return
        for k in range(remaining + 1):
            for combo in combinations(range(i, L + p + 1), k):
                row = tuple(sorted(combo, reverse=True))
                for rest in rows_gen(i + 1, remaining - k):
                    yield [row] + rest
    for rows in rows_gen(1, L):
        word = [x for r in rows for x in r]
        if is_reduced(word):
            T = RCGraph(rows)
            if maxd(T.perm) <= p:
                yield T

tbl = defaultdict(list); n = 0
for p in range(2, 5):
    for L in range(1, 7 if p < 4 else 6):
        for A in rc_graphs(p, L):
            if not is_hw(A, p):
                continue
            n += 1
            rows = [tuple(r) for r in A]
            key = (p, tuple(rows[1:]), len(rows[0]), extwt(A, p))
            tbl[key].append(rows[0])
amb = {k: v for k, v in tbl.items() if len(v) > 1}
print("hw graphs", n, "keys", len(tbl), "ambiguous", len(amb))
for k, v in list(amb.items())[:8]:
    print("  ", k, v)
