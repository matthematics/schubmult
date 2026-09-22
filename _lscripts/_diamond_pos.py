"""Position-based diamond: beta at word position q. Let C = chain positions of L on T.
If q not in C: is L(beta T) == bump at position q of L(T)?  (position-based transport)
If q in C: classify."""
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

def Lchain(word, p):
    w = perm_of(word); m = maxd(w)
    if m <= p:
        return None, None
    q = next(q for q in range(len(word)) if rroot(word, q) == (m, m + 1))
    return ltb_down(word, q, return_chain=True)

def bump_at(word, seq, p, q):
    w = perm_of(word); a, b = rroot(word, q)
    if a <= p or not is_cover(w, a, b):
        return None
    nw = ltb_down(word, q)
    if nw is None or not valid_rc(nw, seq):
        return None
    return nw

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
                LT, C = Lchain(word, p)
                if LT is None:
                    continue
                qL = C[0]
                for q in range(len(word)):
                    if q == qL:
                        continue
                    bT = bump_at(word, seq, p, q)
                    if bT is None:
                        continue
                    nwb, Cb = ltb_down(word, q, return_chain=True)
                    LbT, _ = Lchain(bT, p)
                    inC = q in C
                    if not inC:
                        bLT = bump_at(LT, seq, p, q)
                        ok = (LbT is not None and bLT is not None and LbT == bLT)
                        overlap = bool(set(C) & set(Cb))
                        stats[("q not in C(L)", "L(bT)=bump_q(LT)", ok, "chains overlap", overlap)] += 1
                        if not ok and len(ex) < 6:
                            ex.append(("q notin C", rows[:p], p, q, C, Cb, rroot(word, q)))
                    else:
                        # q in chain of L: L's chain passes through beta's crossing
                        if LbT == LT:
                            stats[("q in C(L)", "L(bT)=L(T)")] += 1
                        else:
                            # is bT = bump at some position of LT?
                            hits = [r for r in range(len(word)) if bump_at(LT, seq, p, r) == bT]
                            hits2 = [r for r in range(len(word)) if LbT is not None and bump_at(LT, seq, p, r) == LbT]
                            stats[("q in C(L)", "bT=bump_r(LT)" if hits else ("L(bT)=bump_r(LT)" if hits2 else "other"))] += 1
                            if not hits and not hits2 and len(ex) < 8:
                                ex.append(("q in C", rows[:p], p, q, C, Cb))
for k, v in sorted(stats.items(), key=lambda x: -x[1]):
    print(v, k)
for e in ex:
    print("  ", e)
