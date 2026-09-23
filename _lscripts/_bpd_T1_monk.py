"""Monk preimage analysis for T1.

R: RC graph, x = leftmost cross of the last nonempty row alpha, at column j, letter a = alpha + j - 1.
R^- = R \\ x in PD(pi), R in PD(pi s_a).  For j >= 2:
  D0 = R^- + (alpha, j-1) reduced  <=>  pi(a-1) < pi(a);  then R = m_{a-1,a}(D0).
Otherwise find the Monk preimage D0' of R (m_{s,a}(D0') = R) by brute force and report.
Also test BPD-side confinement: phi^{-1}(R)_{<alpha} == phi^{-1}(D0')_{<alpha}."""
import sys
from collections import Counter
from itertools import permutations

import numpy as np

from _lb_common import rroots
from schubmult import BPD, Permutation, RCGraph
from schubmult.utils.perm_utils import is_reduced

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5


def cells(rows):
    """(row, col) 1-based cells with letters; rows are tuples of letters (decreasing)."""
    out = []
    for i, r in enumerate(rows, start=1):
        for l in r:
            out.append((i, l - i + 1))
    return out


def rows_from_cells(cs, height):
    rows = [[] for _ in range(height)]
    for i, j in cs:
        rows[i - 1].append(i + j - 1)
    return [tuple(sorted(r, reverse=True)) for r in rows]


def word_of(cs):
    cs = sorted(cs, key=lambda c: (c[0], -c[1]))
    return [i + j - 1 for i, j in cs], cs


def monk_m(cs, s, beta):
    """GH Def 4.2 m_{s,beta} on the pipe dream with cell set cs (reduced, perm = pi t_{s,beta})."""
    cs = set(cs)
    word, order = word_of(cs)
    rts = rroots(word)
    k = [q for q in range(len(word)) if set(rts[q]) == {s, beta}]
    assert len(k) == 1, (cs, s, beta, rts)
    i, j = order[k[0]]
    cs.remove((i, j))
    while True:
        jj = j + 1
        while (i, jj) in cs:
            jj += 1
        cs.add((i, jj))
        word, order = word_of(cs)
        if is_reduced(word):
            return frozenset(cs)
        rts = rroots(word)
        kk = order.index((i, jj))
        partners = [q for q in range(len(word)) if q != kk and set(rts[q]) == set(rts[kk])]
        assert len(partners) == 1
        i, j = order[partners[0]]
        cs.remove((i, j))


def bpd_of(cs, N):
    return BPD.from_rc_graph(RCGraph(rows_from_cells(cs, N))).resize(N)


stats = Counter()
shown = 0
for N0 in range(2, NMAX + 1):
    for a_ in permutations(range(1, N0 + 1)):
        w = Permutation(list(a_))
        if w.inv == 0:
            continue
        for R in RCGraph.all_rc_graphs(w, N0):
            rows = [tuple(r) for r in R]
            cs = cells(rows)
            alpha = max(i for i, _ in cs)
            if alpha == 1:
                continue
            j = min(jj for i, jj in cs if i == alpha)
            if j == 1:
                stats["j=1 (x_alpha)"] += 1
                continue
            a = alpha + j - 1
            Rm = [c for c in cs if c != (alpha, j)]
            pi = Permutation.ref_product(*word_of(Rm)[0]) if Rm else Permutation([])
            N = N0 + 2
            D0 = Rm + [(alpha, j - 1)]
            if is_reduced(word_of(D0)[0]):
                stats["D0 reduced"] += 1
                img = monk_m(D0, a - 1, a)
                stats[("m_{a-1,a}(D0)==R", img == frozenset(cs))] += 1
                # BPD confinement check
                B = bpd_of(cs, N)
                B0 = bpd_of(D0, N)
                ok = np.array_equal(B._grid[: alpha - 1], B0._grid[: alpha - 1])
                stats[("confined (D0 reduced)", ok)] += 1
            else:
                stats["D0 not reduced"] += 1
                # brute-force Monk preimage: s < a with pi t_{s,a} > pi cover
                found = []
                for s in range(1, a):
                    if pi[s - 1] < pi[a - 1] and all(not (pi[s - 1] < pi[k - 1] < pi[a - 1]) for k in range(s + 1, a)):
                        sigma = pi.swap(s - 1, a - 1)
                        wt = tuple(len(r) for r in rows)
                        for D in RCGraph.all_rc_graphs(sigma, N0):
                            if tuple(len(r) for r in D) != wt:
                                continue
                            Dcs = cells([tuple(r) for r in D])
                            if monk_m(Dcs, s, a) == frozenset(cs):
                                found.append((s, frozenset(Dcs)))
                stats[("preimage count", len(found))] += 1
                if found:
                    s, Dcs = found[0]
                    B = bpd_of(cs, N)
                    B0 = bpd_of(Dcs, N)
                    ok = np.array_equal(B._grid[: alpha - 1], B0._grid[: alpha - 1])
                    stats[("confined (D0' case)", ok)] += 1
                    same_top = {c for c in Dcs if c[0] < alpha} == {c for c in cs if c[0] < alpha}
                    stats[("D0' agrees with R above alpha", same_top)] += 1
                    if shown < 4 and not same_top:
                        shown += 1
                        print("R", rows, "alpha", alpha, "j", j, "a", a, "pi", pi, "s", s)
                        print("  D0'", rows_from_cells(Dcs, N0))
for k, v in sorted(stats.items(), key=str):
    print(v, k)
