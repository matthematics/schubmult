"""Verify Lemma `forestlift` of writing/schubert_bialgebra_v2.tex.

Statement: let w be a permutation with maxd(w) <= h, and let B, B' be bounded RC
graphs with h rows, permutation w, and empty last row.  If
    Z(B) =_forest Z(B')
then B =_forest B'.

Equivalently: on RC_h(w)^0, the partition induced by  B |-> P(word(Z(B)))
refines the partition induced by  B |-> P(word(B)).

We also record the converse direction and the finer statement that the fibers of
Z over a single forest class are exactly the forest classes upstairs.
"""

import sys
from collections import defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(perm):
    return len(perm.trimcode)


def rc0(perm, h):
    """Bounded RC graphs with h rows, permutation perm, empty last row."""
    out = []
    for R in RCGraph.all_rc_graphs(perm, h):
        if len(R) != h:
            continue
        if len(R[-1]) != 0:
            continue
        out.append(R)
    return out


def main(N=5):
    total = 0
    fail_lift = 0
    fail_conv = 0
    examples = []
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            d = maxd(w)
            if d == 0:
                continue
            for h in range(d, n + 1):
                graphs = rc0(w, h)
                if not graphs:
                    continue
                # group by downstairs forest invariant, and by upstairs one
                down = defaultdict(set)
                up = defaultdict(set)
                for B in graphs:
                    Z = B.zero_out_last_row()
                    down[Z.forest_invariant].add(B)
                    up[B.forest_invariant].add(B)
                total += len(graphs)
                # forestlift: each down-class is inside a single up-class
                for key, blk in down.items():
                    invs = {B.forest_invariant for B in blk}
                    if len(invs) > 1:
                        fail_lift += 1
                        if len(examples) < 5:
                            examples.append(("lift", w, h, sorted(str(b) for b in blk)))
                # converse: each up-class maps into a single down-class?
                for key, blk in up.items():
                    invs = {B.zero_out_last_row().forest_invariant for B in blk}
                    if len(invs) > 1:
                        fail_conv += 1
    print(f"N={N}: instances={total} forestlift failures={fail_lift} converse failures={fail_conv}")
    for e in examples:
        print(e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
