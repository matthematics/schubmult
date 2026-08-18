"""Test the ROOT-based account of antidiagonal moves and of Z.

User's intuition: think of Z via the insertion algorithm, root-based rather
than word-based.  An antidiagonal move is the deletion of a ROOT from row i and
its reinsertion in row i'.  The root, not the letter, is the invariant.  Since
right roots only affect rows at or below, and i' < i, whether the insertion
negates the root in row i governs what happens in row i' as well.

This predicts exactly what was observed: upstairs the ROWS are preserved while
the LETTER changes, because the letter of a given root differs between R and B.

Tests, all on the NONTRIVIAL branch h = maxd(w):

 (A) downstairs root invariance: for the move (i,j) -> (i',j'), is
        rtt_{R_1}(i,j) == rtt_{R_2}(i',j')  ?
 (B) upstairs, is B_2 obtained from B_1 by deleting the crossing carrying that
     same root and reinserting it in row i'?  We locate the root in B_1 and ask
     for its row.
 (C) the root-level intertwining: is the root moved upstairs the same as the
     root moved downstairs?
"""

import sys
from collections import Counter, defaultdict
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


def antidiag_move(c1, c2):
    gone, came = c1 - c2, c2 - c1
    if len(gone) != 1 or len(came) != 1:
        return None
    (i, j), = gone
    (ip, jp), = came
    if i + j != ip + jp:
        return None
    return (i, j, ip, jp)


def root_map(R):
    """crossing -> root, for every crossing of R"""
    return {(i, j): R.right_root_at(i, j) for (i, j) in crossings(R)}


def main(N=6):
    pairs = 0
    A_ok = 0
    B_found = 0
    B_rowmatch = 0
    C_ok = 0
    rowpat = Counter()
    multi = 0
    multi_rootok = 0
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
                    mv = antidiag_move(c1, crossings(R2))
                    if mv is None:
                        continue
                    i, j, ip, jp = mv
                    pairs += 1

                    alpha1 = R1.right_root_at(i, j)
                    alpha2 = R2.right_root_at(ip, jp)
                    if alpha1 == alpha2:
                        A_ok += 1

                    B1, B2 = Zinv[R1], Zinv[R2]
                    d1, d2 = crossings(B1), crossings(B2)
                    rm1 = root_map(B1)
                    rm2 = root_map(B2)

                    # locate the moved root upstairs
                    src = [c for c, a in rm1.items() if a == alpha1]
                    if len(src) == 1:
                        B_found += 1
                        (bi, bj), = src
                        tgt = [c for c, a in rm2.items() if a == alpha1]
                        if len(tgt) == 1:
                            (ti, tj), = tgt
                            rowpat[(bi - i, ti - ip)] += 1
                            if (bi, ti) == (i, ip):
                                B_rowmatch += 1
                            # root-level move: B2 == B1 with that crossing relocated
                            if (d1 - {(bi, bj)}) | {(ti, tj)} == d2:
                                C_ok += 1
                            elif len(examples) < 5:
                                examples.append(
                                    (tuple(w[:n]), h, f"down root {alpha1} rows {i}->{ip}",
                                     f"up rows {bi}->{ti}",
                                     [tuple(r) for r in R1], [tuple(r) for r in R2],
                                     [tuple(r) for r in B1], [tuple(r) for r in B2])
                                )
                    if len(d1 - d2) > 1:
                        multi += 1
                        if len(src) == 1:
                            (bi, bj), = src
                            tgt = [c for c, a in rm2.items() if a == alpha1]
                            if len(tgt) == 1:
                                (ti, tj), = tgt
                                if (d1 - {(bi, bj)}) | {(ti, tj)} == d2:
                                    multi_rootok += 1

    print(f"N={N}  (h = maxd(w), nontrivial branch)")
    print(f"  same-permutation antidiagonal pairs: {pairs}")
    print()
    print(f"  (A) downstairs root invariance  rtt_R1(i,j) == rtt_R2(i',j'): {A_ok}/{pairs}")
    print(f"  (B) that root occurs exactly once in B_1:                    {B_found}/{pairs}")
    print(f"      ... and its rows upstairs match (i, i'):                 {B_rowmatch}/{pairs}")
    print(f"      (row_up - row_down, tgt_up - tgt_down) histogram: {dict(sorted(rowpat.items()))}")
    print(f"  (C) B_2 == B_1 with that ROOT relocated:                     {C_ok}/{pairs}")
    print()
    print(f"  instances where >1 crossing moves upstairs: {multi}, of which root account still works: {multi_rootok}")
    for e in examples:
        print("\n   ROOT ACCOUNT FAILS:")
        print(f"     w={e[0]} h={e[1]}  {e[2]}   {e[3]}")
        print(f"     R1={e[4]}  R2={e[5]}")
        print(f"     B1={e[6]}  B2={e[7]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
