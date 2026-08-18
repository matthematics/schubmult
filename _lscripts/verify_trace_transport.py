"""Does the traced labelling explain the lifting of antidiagonal moves?

The tracing gives a canonical ROW-PRESERVING bijection

    phi_B :  crossings(B)  ->  crossings(Z(B)).

For an antidiagonal move downstairs, R2 = R1 with the crossing at (i,j) moved to
(i',j'), and B1, B2 the unique preimages, the natural claim is

    B2  =  B1  -  {phi_{B1}^{-1}(i,j)}  +  {phi_{B2}^{-1}(i',j')},

i.e. the move upstairs is the transport of the move downstairs along the
tracing.  Since phi preserves rows, this automatically moves a crossing from row
i to row i', which is exactly the pattern observed (rows preserved, letters not).

We also ask whether phi is just the naive row-order bijection, in which case the
tracing adds nothing.
"""

import sys
from collections import Counter
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph

sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from trace_zero_labels import _cross, maxd, rc0, traced_zero


def row_order_bij(B, R):
    """naive bijection: k-th crossing of row i in B <-> k-th of row i in R"""
    m = {}
    for i in range(1, len(R) + 1):
        cb = [(i, L - i + 1) for L in B[i - 1]]
        cr = [(i, L - i + 1) for L in R[i - 1]]
        if len(cb) != len(cr):
            return None
        for a, b in zip(cb, cr):
            m[b] = a          # crossings(R) -> crossings(B)
    return m


def single_move(c1, c2):
    gone, came = c1 - c2, c2 - c1
    if len(gone) != 1 or len(came) != 1:
        return None
    (i, j), = gone
    (ip, jp), = came
    return (i, j, ip, jp)


def main(N=5):
    same_as_roworder = 0
    total_graphs = 0
    pairs = 0
    transport_ok = 0
    transport_roworder_ok = 0
    src_row_ok = 0
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
            if not graphs:
                continue

            info = {}
            for B in graphs:
                R, labels = traced_zero(B)
                inv = {dst: src for dst, src in labels.items()}   # crossings(R) -> crossings(B)
                info[R] = (B, inv)
                ro = row_order_bij(B, R)
                total_graphs += 1
                if ro is not None and ro == inv:
                    same_as_roworder += 1

            Rs = list(info)
            for R1 in Rs:
                c1 = _cross(R1)
                for R2 in Rs:
                    if R1 is R2 or R1.perm != R2.perm:
                        continue
                    mv = single_move(c1, _cross(R2))
                    if mv is None:
                        continue
                    i, j, ip, jp = mv
                    if i + j != ip + jp:
                        continue
                    pairs += 1
                    B1, inv1 = info[R1]
                    B2, inv2 = info[R2]
                    src = inv1.get((i, j))
                    tgt = inv2.get((ip, jp))
                    if src is None or tgt is None:
                        continue
                    if src[0] == i and tgt[0] == ip:
                        src_row_ok += 1
                    d1, d2 = _cross(B1), _cross(B2)
                    if (d1 - {src}) | {tgt} == d2:
                        transport_ok += 1
                    elif len(examples) < 5:
                        examples.append(
                            (tuple(w[:n]), h, f"down ({i},{j})->({ip},{jp})",
                             f"traced src {src} tgt {tgt}",
                             [tuple(r) for r in R1], [tuple(r) for r in R2],
                             [tuple(r) for r in B1], [tuple(r) for r in B2])
                        )
                    ro1 = row_order_bij(B1, R1)
                    ro2 = row_order_bij(B2, R2)
                    if ro1 and ro2:
                        s2, t2 = ro1.get((i, j)), ro2.get((ip, jp))
                        if s2 and t2 and (d1 - {s2}) | {t2} == d2:
                            transport_roworder_ok += 1

    print(f"N={N}  (h = maxd(w), nontrivial branch)")
    print(f"  traced bijection equals the naive row-order bijection: {same_as_roworder}/{total_graphs}")
    print()
    print(f"  antidiagonal pairs: {pairs}")
    print(f"    src in row i and tgt in row i' (automatic if phi preserves rows): {src_row_ok}/{pairs}")
    print(f"    B2 == B1 with traced src replaced by traced tgt:   {transport_ok}/{pairs}")
    print(f"    same with the naive row-order bijection:           {transport_roworder_ok}/{pairs}")
    for e in examples:
        print("\n   TRANSPORT FAILS:")
        print(f"     w={e[0]} h={e[1]}  {e[2]}   {e[3]}")
        print(f"     R1={e[4]}  R2={e[5]}")
        print(f"     B1={e[6]}  B2={e[7]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
