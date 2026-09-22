"""Specific pairs: Y with rows p+2..m empty, last row m empty, maxd = m, row p+1 nonempty.
T = Y_{<=p}, T' = (Z Y)_{<=p}, matched position-wise. For positions k:
(a) cover status of the p-high root agrees between T and T'? (split: k in row 1 / rows>=2)
(b) chains equal as position sets?
(c) results satisfy D^p(beta T) == D^p(beta' T')?  (expected: confluence)"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, valid_rc, maxd
from collections import Counter

def perm_of(word):
    return Permutation.ref_product(*word) if word else Permutation([])
def arr(w, n=14):
    a = list(w); a += list(range(len(a) + 1, n + 1)); return a[:n]
def is_cover(w, a, b):
    ar = arr(w)
    if ar[a - 1] < ar[b - 1]: return False
    return not any(ar[b - 1] < ar[c - 1] < ar[a - 1] for c in range(a + 1, b))
def is_hw(T, p):
    return all(T.raising_operator(i) is None for i in range(1, p))
def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

stats = Counter(); ex = []
for N in range(2, 9):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0: continue
        m = maxd(w)
        for p in range(2, m - 1):
            for Y in RCGraph.all_rc_graphs(w, m):
                rows = [tuple(r) for r in Y]
                if len(rows[p]) == 0 or any(len(rows[i]) for i in range(p + 1, m)): continue
                hw = is_hw(Y, p)
                ZY = Y.zero_out_last_row()
                T = RCGraph(rows[:p]); T2 = RCGraph([tuple(r) for r in ZY][:p])
                wa, sa = T.as_reduced_compatible(); wb, sb = T2.as_reduced_compatible()
                wa, wb = tuple(wa), tuple(wb); sa = tuple(sa)
                assert tuple(sb) == sa
                u, v = perm_of(wa), perm_of(wb)
                DT, DT2 = Dp(T, p), Dp(T2, p)
                stats[("D^p equal", DT == DT2, "hw", hw)] += 1
                for k in range(len(wa)):
                    ra, rb = rroot(wa, k), rroot(wb, k)
                    ca = ra[0] > p and is_cover(u, *ra); cb = rb[0] > p and is_cover(v, *rb)
                    rowk = "row1" if sa[k] == 1 else "row>=2"
                    stats[("cover agree", ca == cb, rowk, "hw", hw)] += 1
                    if not (ca and cb): continue
                    na, Ca = ltb_down(wa, k, return_chain=True); nb, Cb = ltb_down(wb, k, return_chain=True)
                    if na is None or nb is None or not valid_rc(na, sa) or not valid_rc(nb, sa): continue
                    stats[("chains equal", set(Ca) == set(Cb), rowk, "hw", hw)] += 1
                    stats[("roots equal at k", ra == rb, rowk, "hw", hw)] += 1
                    A = RCGraph.from_reduced_compatible(list(na), list(sa)); B = RCGraph.from_reduced_compatible(list(nb), list(sa))
                    stats[("D^p(bT)==D^p(b'T')", Dp(A, p) == Dp(B, p), rowk, "hw", hw)] += 1
                    if set(Ca) != set(Cb) and rowk == "row>=2" and hw and len(ex) < 4:
                        ex.append((rows[:p + 1], p, k, wa, wb, Ca, Cb))
for k, v in sorted(stats.items(), key=lambda x: str(x[0])):
    print(v, k)
for e in ex: print("  ", e)
