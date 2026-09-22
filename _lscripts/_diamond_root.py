"""Refine the generic diamond L(beta T) = beta'(L T): identify beta' root in terms of beta root (a,b) and the
transition permutation sigma = w_T^{-1} w_{LT} (positions)."""
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
    ar = arr(w, max(b + 1, 14))
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
                m = md
                wa = arr(w); wl = arr(perm_of(LT))
                # sigma as a map on positions: w_LT = w_T * sigma  => sigma = w^{-1} w_LT ; as position permutation:
                # position x of w_LT holds value wl[x]; find position y in w with wa[y]==wl[x]
                sigma = {x + 1: wa.index(wl[x]) + 1 for x in range(len(wl))}
                B = bumps(word, seq, p)
                BL = bumps(LT, seq, p)
                for (a, b), bT in B.items():
                    if bT == LT:
                        continue
                    LbT = L(bT, p)
                    m2 = maxd(perm_of(bT))
                    if LbT is None:
                        continue
                    matches = [r for r, x in BL.items() if x == LbT]
                    if not matches:
                        stats[("no beta'", "maxd(bT)==m", m2 == m)] += 1
                        continue
                    rb = matches[0]
                    cand = tuple(sorted((sigma[a], sigma[b])))
                    stats[("beta' == (a,b)", rb == (a, b), "beta' == sigma(a,b)", rb == cand, "maxd(bT)==m", m2 == m, "a==m", a == m)] += 1
                    if rb != (a, b) and rb != cand and len(ex) < 6:
                        ex.append((rows[:p], p, (a, b), rb, cand, "m", m, "w", "".join(map(str, wa[:7])), "wLT", "".join(map(str, wl[:7]))))
for k, v in sorted(stats.items(), key=lambda x: -x[1]):
    print(v, k)
for e in ex:
    print("  ", e)
