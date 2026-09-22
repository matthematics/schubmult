"""(diamond): for T in RC_{<=p} (maxd>p) and a p-high cover bump beta, find minimal (k,l,j) with
L^k(beta T) in Bumps^{<=j}(L^l T), l>=1 (or beta T normal and equal to L^l T). Report distribution and failures."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, valid_rc, maxd
from collections import Counter

def perm_of(word):
    return Permutation.ref_product(*word) if word else Permutation([])

def is_cover(w, a, b):
    arr = list(w) + list(range(len(w) + 1, b + 2))
    if arr[a - 1] < arr[b - 1]:
        return False
    return not any(arr[b - 1] < arr[c - 1] < arr[a - 1] for c in range(a + 1, b))

def L(word, p):
    w = perm_of(word); m = maxd(w)
    if m <= p:
        return None
    q = next(q for q in range(len(word)) if rroot(word, q) == (m, m + 1))
    return ltb_down(word, q)

def bumps(word, seq, p):
    w = perm_of(word); out = set()
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p and is_cover(w, a, b):
            nw = ltb_down(word, q)
            if nw is not None and valid_rc(nw, seq):
                out.add(nw)
    return out

def bump_closure(word, seq, p, j):
    layers = [{word}]
    for _ in range(j):
        nxt = set()
        for x in layers[-1]:
            nxt |= bumps(x, seq, p)
        layers.append(nxt)
    return layers

stats = Counter(); fails = []
KMAX, LMAX_, JMAX = 3, 3, 3
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
                # L^l T
                Ls = [word]
                for _ in range(LMAX_):
                    nx = L(Ls[-1], p)
                    if nx is None:
                        break
                    Ls.append(nx)
                closures = {l: bump_closure(Ls[l], seq, p, JMAX) for l in range(1, len(Ls))}
                for bT in bumps(word, seq, p):
                    if bT == Ls[1]:
                        stats[("beta==L", 0, 1, 0)] += 1; continue
                    found = None
                    cur = bT
                    for k in range(KMAX + 1):
                        if cur is None:
                            break
                        for l in range(1, len(Ls)):
                            for j in range(JMAX + 1):
                                if cur in closures[l][j]:
                                    found = (k, l, j); break
                            if found: break
                        if found: break
                        cur = L(cur, p)
                    if found:
                        stats[found] += 1
                    else:
                        stats["FAIL"] += 1
                        if len(fails) < 5:
                            fails.append((rows[:p], p, bT))
print(sorted(stats.items(), key=lambda x: -x[1]))
for f in fails:
    print("  ", f)
