"""In the clean class (beta at root (a,m), a<m, same maxd), identify beta' root in LT: candidates (a,m+1), (a,m), pipe transport, position transport."""
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
def bumps(word, seq, p):
    w = perm_of(word); out = {}
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p and is_cover(w, a, b):
            res = ltb_down(word, q, return_chain=True)
            if res[0] is not None and valid_rc(res[0], seq):
                out[q] = res
    return out

stats = Counter(); ex = []
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        m = maxd(w); wa = arr(w)
        for p in range(1, m):
            for T in RCGraph.all_rc_graphs(w, m):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, m)):
                    continue
                word, seq = T.as_reduced_compatible(); word = tuple(word); seq = tuple(seq)
                LT, CL = Lchain(word, p)
                if LT is None:
                    continue
                BL = bumps(LT, seq, p); wl = arr(perm_of(LT))
                for qb, (bT, Cb) in bumps(word, seq, p).items():
                    a, b = rroot(word, qb)
                    if b != m or a == m:
                        continue
                    LbT, _ = Lchain(bT, p)
                    hits = [(q, rroot(LT, q)) for q, (x, _) in BL.items() if x == LbT]
                    if not hits:
                        stats["no beta'"] += 1; continue
                    q2, r2 = hits[0]
                    pipes = tuple(sorted((wl.index(wa[a - 1]) + 1, wl.index(wa[m - 1]) + 1)))
                    stats[("beta' root == (a,m+1)", r2 == (a, m + 1), "== (a,m)", r2 == (a, m), "== pipe transport", r2 == pipes, "same position", q2 == qb, "disjoint", not (set(CL) & set(Cb)))] += 1
                    if r2 != (a, m + 1) and len(ex) < 5:
                        ex.append((rows[:p], p, (a, b), r2, pipes, q2, qb))
for k, v in sorted(stats.items(), key=lambda x: -x[1]):
    print(v, k)
for e in ex:
    print("  ", e)
