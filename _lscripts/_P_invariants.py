"""Invariants of P-tableaux under p-high bumps. Test:
 (i) entries <= p of P unchanged (as a set of cells with values)?
 (ii) restriction of P to entries <= p (shape) unchanged?
 (iii) restriction to entries < p+? etc."""
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

def cells_le(P, t):
    return frozenset((i, j, v) for i, row in enumerate(P) for j, v in enumerate(row) if v <= t)

def shape_le(P, t):
    return tuple(sum(1 for v in row if v <= t) for row in P)

random.seed(0)
stats = Counter()
ex = []
for N in range(3, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 8:
            continue
        words = list(reduced_words(w))
        if len(words) > 20:
            words = random.sample(words, 20)
        seenP = set()
        for word in words:
            P = eg_insert(word)
            if P in seenP:
                continue
            seenP.add(P)
            r = rw(P)
            for p in range(1, maxd(w)):
                for k in range(len(r)):
                    a, b = rroot(r, k)
                    if a <= p:
                        continue
                    nr = ltb_down(r, k)
                    if nr is None:
                        continue
                    P2 = eg_insert(nr)
                    stats[("cells<=p same", cells_le(P, p) == cells_le(P2, p))] += 1
                    stats[("shape<=p same", shape_le(P, p) == shape_le(P2, p))] += 1
                    # which entries decremented: values before
                    dec = [r[i] for i in range(len(r)) if r[i] != nr[i]]
                    stats[("min decremented entry > p", min(dec) > p)] += 1
                    if min(dec) <= p and len(ex) < 4:
                        ex.append((P, p, (a, b), P2))
print(stats)
for e in ex:
    print(e)
