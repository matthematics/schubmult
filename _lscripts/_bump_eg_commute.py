"""Does a Little bump at inversion (a,b) commute with Edelman-Greene insertion?
Test: P(ltb_k w) == P(ltb_{k'} rw(P(w))) where k' is the position of rw(P(w)) with the same right root."""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation
from _lb_common import rroot, ltb_down, reduced_words, maxd
from collections import Counter

def eg_insert(word):
    P = []
    for x in word:
        r = 0
        while True:
            if r == len(P):
                P.append([x]); break
            row = P[r]
            bigger = [k for k, y in enumerate(row) if y > x]
            if not bigger:
                row.append(x); break
            k = bigger[0]; y = row[k]
            if y == x + 1 and k > 0 and row[k - 1] == x:
                x = y
            else:
                row[k] = x; x = y
            r += 1
    return tuple(tuple(r) for r in P)

def rw(P):
    return tuple(x for row in reversed(P) for x in row)

random.seed(0)
stats = Counter()
bad = []
for N in range(3, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 8:
            continue
        words = list(reduced_words(w))
        if len(words) > 25:
            words = random.sample(words, 25)
        for word in words:
            P = eg_insert(word)
            r = rw(P)
            assert Permutation.ref_product(*r) == w
            for k in range(len(word)):
                nw = ltb_down(word, k)
                if nw is None:
                    continue
                root = rroot(word, k)
                kp = [i for i in range(len(r)) if rroot(r, i) == root]
                assert len(kp) == 1
                nr = ltb_down(r, kp[0])
                if nr is None:
                    stats["tableau-bump degenerate"] += 1
                    continue
                P1 = eg_insert(nw)
                P2 = eg_insert(nr)
                stats[("P equal", P1 == P2, "bumped rw is tableau word", nr == rw(P2))] += 1
                if P1 != P2 and len(bad) < 5:
                    bad.append((word, k, root, nw, P1, r, nr, P2))
print(stats)
for b in bad:
    print(b)
