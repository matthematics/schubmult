"""The right root-level question.

(A) established: an antidiagonal move preserves the root, i.e. it deletes the
crossing carrying root alpha in row i and reinserts a crossing carrying the SAME
root alpha in row i'.  So the natural class of moves is:

    ROOT MOVE: delete a crossing at (i,j), insert one at (i',j'), such that the
    root is unchanged.

Downstairs, antidiagonal <=> root move (we test the converse too).

The previous test wrongly looked for the DOWNSTAIRS root alpha inside B_1.  Z
changes the permutation, so alpha need not be a root of B_1 at all.  The correct
question is whether the upstairs move is a root move IN ITS OWN RIGHT, and
whether its rows are (i, i').

Tests, on the nontrivial branch h = maxd(w):
  (1) downstairs: antidiagonal => root preserved (re-confirm), and
      root preserved => antidiagonal (converse)
  (2) upstairs: when a single crossing moves, are the rows exactly (i, i')?
  (3) upstairs: is that move root preserving?
  (4) what is the relation between the upstairs root beta and alpha?
  (5) the multi-crossing cases: what happens there?
"""

import sys
from collections import Counter
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(perm):
    return len(perm.trimcode)


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def crossings(R):
    out = set()
    for i, row in enumerate(R, start=1):
        for L in row:
            out.add((i, L - i + 1))
    return out


def single_move(c1, c2):
    gone, came = c1 - c2, c2 - c1
    if len(gone) != 1 or len(came) != 1:
        return None
    (i, j), = gone
    (ip, jp), = came
    return (i, j, ip, jp)


def main(N=6):
    pairs = 0
    A_ok = 0
    conv_tot = conv_ok = 0
    up_single = 0
    up_rows_ok = 0
    up_rootmove = 0
    up_antidiag = 0
    beta_eq_alpha = Counter()
    multi = 0
    multi_detail = Counter()
    examples = []

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            h = maxd(w)
            if h < 2 or (w, h) in seen:
                continue
            seen.add((w, h))

            graphs = rc0(w, h)
            if len(graphs) < 2:
                continue
            Zinv = {}
            for B in graphs:
                Zinv[B.zero_out_last_row()] = B
            Rs = list(Zinv)

            for R1 in Rs:
                c1 = crossings(R1)
                for R2 in Rs:
                    if R1 is R2 or R1.perm != R2.perm:
                        continue
                    mv = single_move(c1, crossings(R2))
                    if mv is None:
                        continue
                    i, j, ip, jp = mv
                    a1 = R1.right_root_at(i, j)
                    a2 = R2.right_root_at(ip, jp)
                    anti = (i + j == ip + jp)

                    # (1) converse direction: root preserved => antidiagonal?
                    if a1 == a2:
                        conv_tot += 1
                        if anti:
                            conv_ok += 1

                    if not anti:
                        continue
                    pairs += 1
                    if a1 == a2:
                        A_ok += 1

                    B1, B2 = Zinv[R1], Zinv[R2]
                    d1, d2 = crossings(B1), crossings(B2)
                    umv = single_move(d1, d2)
                    if umv is None:
                        multi += 1
                        multi_detail[len(d1 - d2)] += 1
                        continue
                    up_single += 1
                    bi, bj, ti, tj = umv
                    if (bi, ti) == (i, ip):
                        up_rows_ok += 1
                    if bi + bj == ti + tj:
                        up_antidiag += 1
                    b1 = B1.right_root_at(bi, bj)
                    b2 = B2.right_root_at(ti, tj)
                    if b1 == b2:
                        up_rootmove += 1
                        beta_eq_alpha[b1 == a1] += 1
                    elif len(examples) < 5:
                        examples.append(
                            (tuple(w[:n]), h, f"down root {a1} rows {i}->{ip}",
                             f"up roots {b1} -> {b2}, rows {bi}->{ti}",
                             [tuple(r) for r in B1], [tuple(r) for r in B2])
                        )

    print(f"N={N}  (h = maxd(w), nontrivial branch)")
    print(f"  downstairs antidiagonal pairs: {pairs}")
    print()
    print(f"  (1) antidiagonal => root preserved:  {A_ok}/{pairs}")
    print(f"      root preserved => antidiagonal:  {conv_ok}/{conv_tot}")
    print()
    print(f"  (2) upstairs is a single crossing move:  {up_single}/{pairs}")
    print(f"      ... with rows exactly (i, i'):       {up_rows_ok}/{up_single}")
    print(f"      ... and antidiagonal:                {up_antidiag}/{up_single}")
    print(f"  (3) upstairs move is ROOT PRESERVING:    {up_rootmove}/{up_single}")
    print(f"  (4) upstairs root equals downstairs root alpha: {dict(beta_eq_alpha)}")
    print()
    print(f"  (5) multi-crossing upstairs: {multi}, sizes {dict(sorted(multi_detail.items()))}")
    for e in examples:
        print("\n   UPSTAIRS NOT ROOT PRESERVING:")
        print(f"     w={e[0]} h={e[1]}  {e[2]}   {e[3]}")
        print(f"     B1={e[4]}  B2={e[5]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
