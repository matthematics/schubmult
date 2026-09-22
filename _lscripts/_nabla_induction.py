"""Tests for an induction on length:
(C1) nabla D^p(T) == D^p(nabla T)   (nabla = delete first letter / first cross in reading order)
(C2) Delta D^p(T) == D^p(Delta T)   (Delta = delete last letter)
(N1) nabla(beta T) reachable from nabla T by <=2 p-high bumps (or equal), for every p-high bump beta
(N2) same for Delta."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, maxd
from collections import Counter

def rc_from(word, seq):
    return RCGraph.from_reduced_compatible(list(word), list(seq))

def Dp_word(word, seq, p):
    """D^p on the RC graph with given word/rows (rows <= p): returns (word, seq)."""
    T = rc_from(word, seq)
    h = max(p, maxd(T.perm))
    T = T.resize(h)
    while len(T) > p:
        T = T.zero_out_last_row()
    return T.as_reduced_compatible()

def high_bumps(word, seq, p):
    out = set()
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p:
            nw = ltb_down(word, q)
            if nw is not None:
                out.add(nw)
    return out

def reach2(word, seq, p):
    r1 = high_bumps(word, seq, p)
    r2 = set()
    for x in r1:
        r2 |= high_bumps(x, seq, p)
    return {word} | r1 | r2

stats = Counter()
bad = []
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
                D = Dp_word(word, seq, p)
                # nabla
                nw, ns = word[1:], seq[1:]
                Dn = Dp_word(nw, ns, p)
                stats[("C1", (tuple(D[0][1:]), tuple(D[1][1:])) == (tuple(Dn[0]), tuple(Dn[1])))] += 1
                # Delta
                dw, ds = word[:-1], seq[:-1]
                Dd = Dp_word(dw, ds, p)
                stats[("C2", (tuple(D[0][:-1]), tuple(D[1][:-1])) == (tuple(Dd[0]), tuple(Dd[1])))] += 1
                # N1, N2
                Rn = reach2(nw, ns, p)
                Rd = reach2(dw, ds, p)
                for bT in high_bumps(word, seq, p):
                    stats[("N1", bT[1:] in Rn)] += 1
                    stats[("N2", bT[:-1] in Rd)] += 1
                    if bT[1:] not in Rn and len(bad) < 3:
                        bad.append(("N1", rows, p, bT))
                    if bT[:-1] not in Rd and len(bad) < 6:
                        bad.append(("N2", rows, p, bT))
print(stats)
for b in bad:
    print(b)
