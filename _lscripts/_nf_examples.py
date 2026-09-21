"""Print examples of P -> NF_p(P) to guess a direct rule."""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation
from _lb_common import rroot, ltb_down, reduced_words, maxd

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

def high_bumps(word, p):
    out = set()
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p:
            nw = ltb_down(word, q)
            if nw is not None:
                out.add(nw)
    return out

def NF(word, p):
    while True:
        nb = high_bumps(word, p)
        if not nb:
            return word
        word = next(iter(nb))

def perm_arr(word, n=8):
    w = Permutation.ref_product(*word) if word else Permutation([])
    arr = list(w); arr += list(range(len(arr) + 1, n + 1))
    return arr[:n]

random.seed(3)
seen = set()
count = 0
for N in range(4, 6):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv < 3 or w.inv > 5:
            continue
        for word in reduced_words(w):
            P = eg_insert(word)
            for p in range(1, maxd(w)):
                if (P, p) in seen:
                    continue
                seen.add((P, p))
                nf = NF(rw(P), p)
                PN = eg_insert(nf)
                if maxd(Permutation.ref_product(*nf)) > p:
                    continue
                if random.random() < 0.08 and count < 25:
                    count += 1
                    print(f"p={p} w={perm_arr(word)} P={P}  ->  wN={perm_arr(nf)} PN={PN}")
