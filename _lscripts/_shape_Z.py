"""Local lemma at the level of Z (one full level of last-descent bumps, until maxd < m):
for beta on T (maxd(w_T)=m>p), classify: Z(bT) vs Z(T): equal; Z(bT)=b'(ZT); bT = b'(ZT); Z^2(bT) ...; also
whether maxd(bT)<m (then Z(bT) := bT)."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, valid_rc, maxd
from collections import Counter, defaultdict

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

def Z(word, p):
    """bumps at last descent until maxd < maxd(word) (or <= p)."""
    m = maxd(perm_of(word))
    if m <= p:
        return None
    cur = word
    while maxd(perm_of(cur)) >= m:
        cur = L(cur, p)
        if cur is None:
            return None
    return cur

def bumps(word, seq, p):
    w = perm_of(word); out = {}
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p and is_cover(w, a, b):
            nw = ltb_down(word, q)
            if nw is not None and valid_rc(nw, seq):
                out[q] = nw
    return out

stats = Counter(); ex = defaultdict(list)
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w); m = md
        for p in range(1, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible(); word = tuple(word); seq = tuple(seq)
                ZT = Z(word, p)
                if ZT is None:
                    continue
                BZ = set(bumps(ZT, seq, p).values())
                BZ2 = set()
                for x in BZ:
                    BZ2 |= set(bumps(x, seq, p).values())
                for qb, bT in bumps(word, seq, p).items():
                    a, b = rroot(word, qb)
                    lowered = maxd(perm_of(bT)) < m
                    ZbT = bT if lowered else Z(bT, p)
                    if ZbT == ZT:
                        s = "Z(bT)=Z(T)"
                    elif ZbT in BZ:
                        s = "Z(bT)=b'(ZT)"
                    elif ZbT in BZ2:
                        s = "Z(bT)=b'b''(ZT)"
                    else:
                        Z2 = Z(ZbT, p)
                        if Z2 is not None and (Z2 == ZT or Z2 in BZ or Z2 in BZ2):
                            s = "Z^2(bT) joins"
                        else:
                            s = "other"
                    key = (s, "beta lowers maxd" if lowered else "same maxd", "a=m" if a == m else "a<m")
                    stats[key] += 1
                    if len(ex[key]) < 2:
                        ex[key].append((rows[:p], p, (a, b)))
for k, v in sorted(stats.items(), key=lambda x: -x[1]):
    print(v, k, ex[k])
