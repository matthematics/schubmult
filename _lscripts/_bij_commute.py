"""Does the bijection of inversion sets induced by a Little bump commute with Little bumps?
Two natural bijections I(w_T) -> I(w_{bT}) for a bump b at position q_b:
  (pos) rt_T(q) -> rt_{bT}(q)   (same word position)
  (pipe) (x,y) -> positions in w_{bT} of the values w_T(x), w_T(y)
Test, for two p-high cover bumps b (pos q_b) and g (pos q_g) of T with q_g not in C_b and q_b not in C_g:
  pos:  bump at position q_g of bT  ==  bump at position q_b of gT ?
  pipe: bump at pipe-transport of rt(q_g) in bT == bump at pipe-transport of rt(q_b) in gT ?"""
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
def bumps(word, seq, p):
    w = perm_of(word); out = {}
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p and is_cover(w, a, b):
            res = ltb_down(word, q, return_chain=True)
            if res[0] is not None and valid_rc(res[0], seq):
                out[q] = res
    return out
def bump_at_root(word, seq, p, root):
    w = perm_of(word)
    for q in range(len(word)):
        if rroot(word, q) == root:
            a, b = root
            if a <= p or not is_cover(w, a, b):
                return None
            nw = ltb_down(word, q)
            return nw if nw is not None and valid_rc(nw, seq) else None
    return None
def pipe_transport(wf, wt, root):
    af, at = arr(wf), arr(wt)
    x, y = af[root[0] - 1], af[root[1] - 1]
    return tuple(sorted((at.index(x) + 1, at.index(y) + 1)))

st = Counter(); ex = {}
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
                B = bumps(word, seq, p)
                qs = sorted(B)
                for i in range(len(qs)):
                    for j in range(i + 1, len(qs)):
                        qb, qg = qs[i], qs[j]
                        bT, Cb = B[qb]; gT, Cg = B[qg]
                        if qg in Cb or qb in Cg:
                            st["starts on other chain (skip)"] += 1; continue
                        disjoint = not (set(Cb) & set(Cg))
                        # pos transport
                        Bb = bumps(bT, seq, p); Bg = bumps(gT, seq, p)
                        pos_ok = (qg in Bb and qb in Bg and Bb[qg][0] == Bg[qb][0])
                        st[("pos", disjoint, pos_ok)] += 1
                        # pipe transport
                        rg = pipe_transport(w, perm_of(bT), rroot(word, qg)); rb = pipe_transport(w, perm_of(gT), rroot(word, qb))
                        x = bump_at_root(bT, seq, p, rg); y = bump_at_root(gT, seq, p, rb)
                        pipe_ok = (x is not None and y is not None and x == y)
                        st[("pipe", disjoint, pipe_ok)] += 1
                        for name, ok in (("pos", pos_ok), ("pipe", pipe_ok)):
                            if not ok and (name, disjoint) not in ex:
                                ex[(name, disjoint)] = (rows[:p], p, rroot(word, qb), rroot(word, qg))
for k, v in sorted(st.items(), key=str):
    print(v, k)
for k, v in ex.items():
    print("counterexample", k, v)
