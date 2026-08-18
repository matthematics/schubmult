"""What IS true upstairs?  Classify the relation between B_1 and B_2 precisely.

Lemma zeroantidiag as stated (B_2 = sigma(B_1), the identical move) is FALSE on
the nontrivial branch.  The data suggest the correct statement preserves the
ROWS but not the letter:

    B_2 is obtained from B_1 by moving a single crossing from row i to row i',
    possibly carrying a different letter than the downstairs move.

We test that, and also re-test the forest bridge on the nontrivial branch,
since its earlier verification was vacuous for the same reason.
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


def single_move(c1, c2):
    """(i, j, i', j') if exactly one crossing moves; else None."""
    gone, came = c1 - c2, c2 - c1
    if len(gone) != 1 or len(came) != 1:
        return None
    (i, j), = gone
    (ip, jp), = came
    return (i, j, ip, jp)


def antidiag_sig(c1, c2):
    mv = single_move(c1, c2)
    if mv is None:
        return None
    i, j, ip, jp = mv
    if i + j != ip + jp:
        return None
    return (i + j - 1, i, ip)


def main(N=6):
    pairs = 0
    up_single = 0
    up_rows_match = 0
    up_antidiag = 0
    up_identical = 0
    ncross = Counter()
    lshift = Counter()
    bridge = Counter()
    examples = []

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            h = maxd(w)
            if h < 2:
                continue
            if (w, h) in seen:
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
                    sig = antidiag_sig(c1, crossings(R2))
                    if sig is None:
                        continue
                    L, i, ip = sig
                    pairs += 1

                    B1, B2 = Zinv[R1], Zinv[R2]
                    d1, d2 = crossings(B1), crossings(B2)
                    ncross[len(d1 - d2)] += 1

                    mv = single_move(d1, d2)
                    if mv is not None:
                        up_single += 1
                        iu, ju, ipu, jpu = mv
                        Lu = iu + ju - 1
                        if (iu, ipu) == (i, ip):
                            up_rows_match += 1
                            lshift[Lu - L] += 1
                        if iu + ju == ipu + jpu:
                            up_antidiag += 1
                            if (Lu, iu, ipu) == sig:
                                up_identical += 1
                        if (iu, ipu) != (i, ip) and len(examples) < 5:
                            examples.append(
                                (tuple(w[:n]), h, f"down rows {i}->{ip} letter {L}",
                                 f"up rows {iu}->{ipu} letter {Lu}",
                                 [tuple(r) for r in R1], [tuple(r) for r in R2],
                                 [tuple(r) for r in B1], [tuple(r) for r in B2])
                            )

                    # forest bridge, on the nontrivial branch
                    down = R1.forest_invariant == R2.forest_invariant
                    up = B1.forest_invariant == B2.forest_invariant
                    bridge[(down, up)] += 1

    print(f"N={N}  (h = maxd(w), NONTRIVIAL branch)")
    print(f"  same-permutation antidiagonal pairs downstairs: {pairs}")
    print(f"  crossings moved upstairs, histogram: {dict(sorted(ncross.items()))}")
    print()
    print(f"  upstairs is a SINGLE crossing move:      {up_single}/{pairs}")
    print(f"  ... with the SAME source/target rows:    {up_rows_match}/{pairs}")
    print(f"  ... and letter shift L_up - L_down:      {dict(sorted(lshift.items()))}")
    print(f"  upstairs move is antidiagonal:           {up_antidiag}/{pairs}")
    print(f"  upstairs move is IDENTICAL to sigma:     {up_identical}/{pairs}")
    print()
    print(f"  FOREST BRIDGE (down preserves, up preserves) -> count: {dict(bridge)}")
    dt = bridge[(True, True)] + bridge[(True, False)]
    ut = bridge[(True, True)] + bridge[(False, True)]
    if dt:
        print(f"     down => up : {bridge[(True, True)]}/{dt}  (failures {bridge[(True, False)]})")
    if ut:
        print(f"     up => down : {bridge[(True, True)]}/{ut}  (failures {bridge[(False, True)]})")
    for e in examples:
        print("\n   ROWS DIFFER:")
        print(f"     w={e[0]} h={e[1]}  {e[2]}   {e[3]}")
        print(f"     R1={e[4]}  R2={e[5]}")
        print(f"     B1={e[6]}  B2={e[7]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
