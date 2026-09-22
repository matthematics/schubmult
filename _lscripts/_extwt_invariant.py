"""Crystal invariant: for hw T in RC_{<=p}, extwt(T) := lex-min weight in the A_{p-1}-component of T
(closure under f_i, i<p). Since Z is a crystal iso, extwt is shared by D^p(Y) and D^p(ZY).
Test: does (p, wt, extwt) determine D^p(T) among hw T? And (p, w, wt, extwt)?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd
from collections import defaultdict

def is_hw(T, p):
    return all(T.raising_operator(i) is None for i in range(1, p))

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

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
    return best, len(seen)

t1 = defaultdict(set); t2 = defaultdict(set)
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(2, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)) or not is_hw(T, p):
                    continue
                wt = tuple(len(r) for r in rows[:p])
                e, sz = extwt(T, p)
                D = Dp(T, p)
                eD, szD = extwt(D, p)
                assert e == eD and sz == szD, (rows, p, e, eD)
                t1[(p, wt, e)].add(tuple(tuple(r) for r in D))
                t2[(p, tuple(w), wt, e)].add(tuple(tuple(r) for r in D))
for name, tb in (("(p,wt,extwt)", t1), ("(p,w,wt,extwt)", t2)):
    amb = {k: v for k, v in tb.items() if len(v) > 1}
    print(name, "keys", len(tb), "ambiguous", len(amb))
    for k, v in list(amb.items())[:4]:
        print("   ", k, v)
