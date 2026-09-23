"""(Conf-L) For every RC graph X (rows <= n, padded) and every cross c of X with right root (s, beta), s < beta,
let X' = leftward Little bump of X at c (Monk preimage).  Claim: phi^{-1}(X) and phi^{-1}(X') agree in rows
< tau(beta) (non-degenerate: no cross removed) resp. rows < s (degenerate: cross pushed off column 1),
where tau = perm(X).  Also record the sharper bound rows <= tau(beta)-1 vs actual first differing row."""
import sys
from collections import Counter
from itertools import permutations

import numpy as np

from _lb_common import rroots
from schubmult import BPD, Permutation, RCGraph
from schubmult.utils.perm_utils import is_reduced

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5


def cells(rows):
    return [(i, l - i + 1) for i, r in enumerate(rows, start=1) for l in r]


def rows_from_cells(cs, height):
    rows = [[] for _ in range(height)]
    for i, j in cs:
        rows[i - 1].append(i + j - 1)
    return [tuple(sorted(r, reverse=True)) for r in rows]


def word_of(cs):
    cs = sorted(cs, key=lambda c: (c[0], -c[1]))
    return [i + j - 1 for i, j in cs], cs


def left_bump(cs, cell):
    """Inverse Monk move: remove cell, insert at nearest elbow to the left in its row, cascade.
    Returns (new cell set, degenerate flag)."""
    cs = set(cs)
    i, j = cell
    while True:
        cs.remove((i, j))
        jj = j - 1
        while jj >= 1 and (i, jj) in cs:
            jj -= 1
        if jj == 0:
            return frozenset(cs), True
        cs.add((i, jj))
        word, order = word_of(cs)
        if is_reduced(word):
            return frozenset(cs), False
        rts = rroots(word)
        kk = order.index((i, jj))
        partners = [q for q in range(len(word)) if q != kk and set(rts[q]) == set(rts[kk])]
        if len(partners) != 1:
            # fall back: the unique other deletable position
            partners = [q for q in range(len(word)) if q != kk and is_reduced(word[:q] + word[q + 1 :])]
        assert len(partners) == 1, (word, kk, rts)
        i, j = order[partners[0]]


def bpd_of(cs, N):
    return BPD.from_rc_graph(RCGraph(rows_from_cells(cs, N))).resize(N)


stats = Counter()
ex = None
for N0 in range(2, NMAX + 1):
    for a_ in permutations(range(1, N0 + 1)):
        w = Permutation(list(a_))
        if w.inv == 0:
            continue
        for X in RCGraph.all_rc_graphs(w, N0):
            cs = cells([tuple(r) for r in X])
            word, order = word_of(cs)
            rts = rroots(word)
            N = N0 + 2
            B = bpd_of(cs, N)
            for k, c in enumerate(order):
                s, beta = sorted(rts[k])
                # only covers: deleting c must leave a reduced word
                if not is_reduced(word[:k] + word[k + 1 :]):
                    stats["skipped (not a cover)"] += 1
                    continue
                Xp, degen = left_bump(cs, c)
                bound = s
                Bp = bpd_of(Xp, N)
                ok = np.array_equal(B._grid[: bound - 1], Bp._grid[: bound - 1])
                stats[("degenerate" if degen else "nondegenerate", ok)] += 1
                if not ok and ex is None:
                    ex = (rows_from_cells(cs, N0), c, (s, beta), bound, rows_from_cells(Xp, N0))
                # how tight: first differing row
                diff = [r for r in range(N) if not np.array_equal(B._grid[r], Bp._grid[r])]
                first = diff[0] + 1 if diff else None
                if first is not None:
                    stats[("slack first_diff - bound", first - bound)] += 1
for k, v in sorted(stats.items(), key=str):
    print(v, k)
print("counterexample:", ex)
