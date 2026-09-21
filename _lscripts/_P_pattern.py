"""Is the order pattern of the cells of P (sign of differences) preserved by p-high bumps on P?"""
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

def sgn(x):
    return (x > 0) - (x < 0)

def pat(P):
    cells = [v for row in P for v in row]
    n = len(cells)
    return tuple(sgn(cells[i] - cells[j]) for i in range(n) for j in range(i + 1, n))

random.seed(0)
stats = Counter(); ex = []
seen = set()
for N in range(3, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 8:
            continue
        words = list(reduced_words(w))
        if len(words) > 10:
            words = random.sample(words, 10)
        for word in words:
            P = eg_insert(word)
            if P in seen:
                continue
            seen.add(P)
            r = rw(P)
            assert eg_insert(r) == P
            for q in range(len(r)):
                a, b = rroot(r, q)
                nr = ltb_down(r, q)
                if nr is None:
                    continue
                P2 = eg_insert(nr)
                same = pat(P2) == pat(P)
                for p in range(1, maxd(w)):
                    stats[("p-high" if a > p else "low", same)] += 1
                if not same and a > 2 and len(ex) < 5:
                    ex.append((P, (a, b), P2))
print(stats)
for e in ex:
    print(e)
