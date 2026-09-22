"""Monotonicity of low-low inversions under p-high bumps and along D^p."""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, reduced_words, maxd
from collections import Counter

def perm_arr(word, n=12):
    w = Permutation.ref_product(*word) if word else Permutation([])
    a = list(w); a += list(range(len(a) + 1, n + 1)); return a[:n]

def inv_ll(word, p):
    a = perm_arr(word)
    return frozenset((i, j) for i in range(p) for j in range(i + 1, p) if a[i] > a[j])

random.seed(0)
stats = Counter(); ex = []
for N in range(3, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 8:
            continue
        words = list(reduced_words(w))
        if len(words) > 12:
            words = random.sample(words, 12)
        for word in words:
            for p in range(2, maxd(w)):
                I0 = inv_ll(word, p)
                for q in range(len(word)):
                    a, b = rroot(word, q)
                    if a <= p:
                        continue
                    nw = ltb_down(word, q)
                    if nw is None:
                        continue
                    I1 = inv_ll(nw, p)
                    stats[("subset", I0 <= I1, "superset", I0 >= I1)] += 1
                    if not (I0 <= I1) and len(ex) < 4:
                        ex.append((word, p, (a, b), nw, perm_arr(word, 7), perm_arr(nw, 7)))
print(stats)
for e in ex:
    print(e)

# along D^p on RC graphs
def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X
st2 = Counter()
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(2, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word = tuple(T.as_reduced_compatible()[0])
                D = tuple(Dp(T, p).as_reduced_compatible()[0])
                st2[("subset", inv_ll(word, p) <= inv_ll(D, p))] += 1
print("along D^p:", st2)
