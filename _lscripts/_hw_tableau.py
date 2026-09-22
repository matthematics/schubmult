"""Highest-weight T in RC_{<=p}: (H1) word(T) equals a canonical reading word of P(T)?
(H2) tabulate P(T) -> P(D^p T)."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd
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

def colword(P):
    ncols = len(P[0]) if P else 0
    return tuple(P[i][j] for j in range(ncols - 1, -1, -1) for i in range(len(P)) if j < len(P[i]))

def rowword(P):
    return tuple(x for row in reversed(P) for x in row)

def is_hw(T, p):
    return all(T.raising_operator(i) is None for i in range(1, p))

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

def arr(w, n=8):
    a = list(w); a += list(range(len(a) + 1, n + 1)); return tuple(a[:n])

h1 = Counter(); rows_out = []
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(2, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)) or not is_hw(T, p):
                    continue
                word = tuple(T.as_reduced_compatible()[0])
                # EG conventions: insert left-to-right or reversed
                for conv, wd in (("L2R", word), ("rev", tuple(reversed(word)))):
                    P = eg_insert(wd)
                    h1[(conv, "word==colword", wd == colword(P))] += 1
                    h1[(conv, "word==rowword", wd == rowword(P))] += 1
                    h1[(conv, "shape==wt", tuple(len(r) for r in P) == tuple(len(r) for r in rows[:p] if len(r)))] += 1
                D = Dp(T, p)
                if w.inv <= 5 and p == 2:
                    PT = eg_insert(tuple(reversed(word)))
                    PD = eg_insert(tuple(reversed(tuple(D.as_reduced_compatible()[0]))))
                    rows_out.append((rows[:p], arr(w), PT, [tuple(r) for r in D], arr(D.perm), PD))
print(h1)
rows_out.sort(key=lambda r: (sum(len(x) for x in r[0]), r[0]))
for r in rows_out[:60]:
    print(f"T={r[0]} w={r[1]} P={r[2]} -> D={r[3]} w'={r[4]} P'={r[5]}")
