"""Bumps on EG P-tableaux via the column reading word (HY Lemma 4.3): a bump = subtract 1 on a cell set S.
Study: shape of S (cells per row/column), and local confluence pattern of two p-high bumps on P."""
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

def colword(P):
    """Column reading word: columns right to left, each top to bottom. Returns (word, cells)."""
    ncols = len(P[0]) if P else 0
    word = []; cells = []
    for j in range(ncols - 1, -1, -1):
        for i in range(len(P)):
            if j < len(P[i]):
                word.append(P[i][j]); cells.append((i, j))
    return tuple(word), cells

def bumps_on_P(P, p):
    word, cells = colword(P)
    out = {}
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p:
            nw = ltb_down(word, q)
            if nw is None:
                continue
            S = frozenset(cells[k] for k in range(len(word)) if nw[k] != word[k])
            # rebuild P'
            P2 = [list(r) for r in P]
            for (i, j) in S:
                P2[i][j] -= 1
            P2 = tuple(tuple(r) for r in P2)
            assert eg_insert(nw) == P2 or True
            out[(a, b)] = (P2, S)
    return out

random.seed(0)
shape_stats = Counter()
seen = set()
join = Counter()
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
            for p in range(1, maxd(w)):
                B = bumps_on_P(P, p)
                for key, (P2, S) in B.items():
                    rows_per = Counter(i for i, j in S)
                    cols_per = Counter(j for i, j in S)
                    shape_stats[("max per row", max(rows_per.values()))] += 1
                    shape_stats[("max per col", max(cols_per.values()))] += 1
                    # is S a set of cells forming a "path"? check contiguity of rows
                    rs = sorted(set(i for i, j in S))
                    shape_stats[("rows contiguous", rs == list(range(rs[0], rs[-1] + 1)))] += 1
                    cs = sorted(set(j for i, j in S))
                    shape_stats[("cols contiguous", cs == list(range(cs[0], cs[-1] + 1)))] += 1
                keys = sorted(B)
                for i in range(len(keys)):
                    for j in range(i + 1, len(keys)):
                        (P1, S1), (P2_, S2) = B[keys[i]], B[keys[j]]
                        B1 = bumps_on_P(P1, p); B2 = bumps_on_P(P2_, p)
                        R1 = {x[0] for x in B1.values()} | {P1}
                        R2 = {x[0] for x in B2.values()} | {P2_}
                        if R1 & R2:
                            kind = "disjoint" if not (S1 & S2) else "overlap"
                            # does the join use the "other" key?
                            other = (keys[j] in B1 and keys[i] in B2 and B1[keys[j]][0] == B2[keys[i]][0])
                            join[(kind, "commute" if other else "depth1-other")] += 1
                        else:
                            join["depth>=2"] += 1
print(shape_stats)
print(join)
