"""Systematic table of P -> NF_p(P) for small P (shape sizes <= 4), to guess a rule."""
import sys
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

def perm_arr(word, n=7):
    w = Permutation.ref_product(*word) if word else Permutation([])
    arr = list(w); arr += list(range(len(arr) + 1, n + 1))
    return tuple(arr[:n])

rows = []
seen = set()
for N in range(2, 7):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 4:
            continue
        for word in reduced_words(w):
            P = eg_insert(word)
            if P in seen:
                continue
            seen.add(P)
            for p in range(1, maxd(w)):
                nf = NF(rw(P), p)
                if maxd(Permutation.ref_product(*nf)) > p:
                    continue
                PN = eg_insert(nf)
                rows.append((p, P, perm_arr(word), PN, perm_arr(nf)))
rows.sort(key=lambda r: (r[0], len(rw(r[1])), r[1]))
for r in rows[:120]:
    print(f"p={r[0]} P={r[1]} w={r[2]} -> PN={r[3]} wN={r[4]}")
