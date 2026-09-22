"""Generic diamond with pipe transport: beta at positions (a,b) of w has pipes (w(a),w(b)).
beta' := bump of LT at the positions holding the same pipes. Classify by pipe-sharing with L's pipes (w(m),w(m+1))."""
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
    if ar[a - 1] < ar[b - 1]:
        return False
    return not any(ar[b - 1] < ar[c - 1] < ar[a - 1] for c in range(a + 1, b))

def L(word, p):
    w = perm_of(word); m = maxd(w)
    if m <= p:
        return None
    q = next(q for q in range(len(word)) if rroot(word, q) == (m, m + 1))
    return ltb_down(word, q)

def bumps(word, seq, p):
    w = perm_of(word); out = {}
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p and is_cover(w, a, b):
            nw = ltb_down(word, q)
            if nw is not None and valid_rc(nw, seq):
                out[(a, b)] = nw
    return out

def transport(word_from, word_to, a, b):
    wf = arr(perm_of(word_from)); wt = arr(perm_of(word_to))
    x, y = wf[a - 1], wf[b - 1]
    pa, pb = wt.index(x) + 1, wt.index(y) + 1
    return tuple(sorted((pa, pb)))

stats = Counter(); ex = []
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(1, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible(); word = tuple(word); seq = tuple(seq)
                LT = L(word, p)
                if LT is None:
                    continue
                m = md; wa = arr(w)
                Lpipes = {wa[m - 1], wa[m]}
                B = bumps(word, seq, p); BL = bumps(LT, seq, p)
                for (a, b), bT in B.items():
                    if bT == LT:
                        continue
                    share = len({wa[a - 1], wa[b - 1]} & Lpipes)
                    tr = transport(word, LT, a, b)
                    LbT = L(bT, p)
                    # candidate joins
                    if LbT is not None and tr in BL and BL[tr] == LbT:
                        kind = "L(bT)=beta'(LT), beta'=transport"
                    elif LbT is not None and LbT == LT:
                        kind = "L(bT)=L(T)"
                    elif bT in BL.values():
                        r = [r for r, x in BL.items() if x == bT][0]
                        kind = "bT=beta'(LT)" + (" transport" if r == tr else " other")
                    elif LbT is not None and LbT in BL.values():
                        kind = "L(bT)=beta'(LT), beta' other root"
                    else:
                        L2 = L(LbT, p) if LbT is not None else None
                        if L2 is not None and (L2 == LT or L2 in BL.values()):
                            kind = "L^2(bT) in {LT} u beta'(LT)"
                        else:
                            kind = "other"
                    stats[(kind, "shared pipes", share, "a==m", a == m)] += 1
                    if kind == "other" and len(ex) < 4:
                        ex.append((rows[:p], p, (a, b)))
for k, v in sorted(stats.items(), key=lambda x: -x[1]):
    print(v, k)
for e in ex:
    print("  ", e)
