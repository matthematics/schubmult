"""Determine the correct Edelman-Greene convention for RC graphs:
T(R) := recording tableau with entry k replaced by the row of the k-th inserted letter (insertion right to left),
optionally with rows relabeled r -> m+1-r.  Check semistandardness and which tableau crystal operator matches
the RC-graph operator e_i of Definition crystal (bracketing on column indices)."""
import sys
from collections import Counter
from itertools import permutations

from schubmult import Permutation, RCGraph

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
M = int(sys.argv[2]) if len(sys.argv) > 2 else 3  # number of rows


def eg_insert(P, x):
    """Insert x into row 0 of P (list of lists, English). Returns position (row, col) of the new box."""
    P = [list(r) for r in P]
    row = 0
    while True:
        if row == len(P):
            P.append([x])
            return P, (row, 0)
        R = P[row]
        if not R or x >= R[-1]:
            R.append(x)
            return P, (row, len(R) - 1)
        j = next(k for k, v in enumerate(R) if x < v)
        if R[j] == x + 1 and j > 0 and R[j - 1] == x:
            x = x + 1  # special bump: row unchanged
        else:
            R[j], x = x, R[j]
        row += 1


def eg(word, rows):
    P, Q = [], []
    m = len(word)
    for k in range(m):
        letter = word[m - 1 - k]
        r = rows[m - 1 - k]
        P, (i, j) = eg_insert(P, letter)
        while len(Q) <= i:
            Q.append([])
        assert j == len(Q[i])
        Q[i].append(r)
    return P, Q


def is_ssyt(T):
    for r in T:
        if any(r[k] > r[k + 1] for k in range(len(r) - 1)):
            return False
    for i in range(len(T) - 1):
        for j in range(len(T[i + 1])):
            if T[i][j] >= T[i + 1][j]:
                return False
    return True


def tab_reading_word(T):
    # English: rows bottom to top, left to right
    return [x for r in reversed(T) for x in r]


def tab_op(T, i, kind):
    """Tableau crystal operator e_i / f_i on SSYT T via bracketing on the reading word."""
    cells = [(ri, ci) for ri in reversed(range(len(T))) for ci in range(len(T[ri]))]
    word = [T[ri][ci] for ri, ci in cells]
    # pair i+1 (open) with later i (close)? Standard: scan left to right, i+1 = '(' , i = ')'  for f_i? Use:
    # treat i as ')' and i+1 as '(': unmatched ')' i's and unmatched '(' i+1's.
    stack = []
    unpaired_i = []
    unpaired_ip1 = []
    for pos, x in enumerate(word):
        if x == i + 1:
            stack.append(pos)
        elif x == i:
            if stack:
                stack.pop()
            else:
                unpaired_i.append(pos)
    unpaired_ip1 = stack
    T2 = [list(r) for r in T]
    if kind == "f":
        if not unpaired_i:
            return None
        pos = unpaired_i[-1]  # rightmost unpaired i
        ri, ci = cells[pos]
        T2[ri][ci] = i + 1
    else:
        if not unpaired_ip1:
            return None
        pos = unpaired_ip1[0]  # leftmost unpaired i+1
        ri, ci = cells[pos]
        T2[ri][ci] = i
    return T2


def rc_e(R, i):
    """e_i of Definition crystal on RC graph rows (tuples of letters); i is 1-based row index."""
    rows = [list(r) for r in R]
    if i + 1 > len(rows):
        return None
    # column indices: letter l in row r sits at column l - r + 1
    top = sorted([l - i + 1 for l in rows[i - 1]], reverse=True)  # row i, decreasing column
    bot = set(l - (i + 1) + 1 for l in rows[i])  # row i+1 columns
    unpaired_bot = set(bot)
    for b in top:
        cands = [a for a in unpaired_bot if a > b]
        if cands:
            unpaired_bot.remove(min(cands))
    if not unpaired_bot:
        return None
    a = max(unpaired_bot)
    topcols = set(top)
    ap = a + 1
    while ap in topcols:
        ap += 1
    new_rows = [list(r) for r in rows]
    new_rows[i].remove(a + (i + 1) - 1)
    new_rows[i - 1].append(ap + i - 1)
    new_rows[i - 1].sort(reverse=True)
    return [tuple(r) for r in new_rows]


def data(rows_letters, relabel):
    rows = [tuple(r) for r in rows_letters]
    word = [l for r in rows for l in r]
    rowidx = [ri + 1 for ri, r in enumerate(rows) for _ in r]
    if relabel:
        rowidx = [M + 1 - r for r in rowidx]
    return eg(word, rowidx)


stats = Counter()
for N0 in range(2, NMAX + 1):
    for a in permutations(range(1, N0 + 1)):
        w = Permutation(list(a))
        if w.inv == 0:
            continue
        for R in RCGraph.all_rc_graphs(w, M):
            rows = [tuple(r) for r in R]
            for relabel in (False, True):
                P, T = data(rows, relabel)
                stats[("ssyt", relabel, is_ssyt(T))] += 1
                for i in range(1, M):
                    Re_ = RCGraph(rows).raising_operator(i); Re = None if Re_ is None else [tuple(r) for r in Re_]
                    if Re is None:
                        continue
                    if not RCGraph(Re).is_valid if hasattr(RCGraph(Re), "is_valid") else False:
                        stats["invalid e_i output"] += 1
                    P2, T2 = data(Re, relabel)
                    stats[("P preserved", relabel, P2 == P)] += 1
                    for kind in ("e", "f"):
                        for idx in (i, M - i):
                            if idx < 1:
                                continue
                            Tt = tab_op(T, idx, kind)
                            stats[("match", relabel, kind, "i" if idx == i else "M-i", Tt == T2)] += 1
for k, v in sorted(stats.items(), key=str):
    print(v, k)
