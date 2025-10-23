"""
Check behaviour of RCGraph product vs. Edelman–Greene Q-tableau.

Usage:
    python -m schubmult.scripts.check_q_tableau N

This script is exploratory / best-effort.
"""
import sys
import itertools
from collections import Counter
from itertools import zip_longest

from schubmult.perm_lib import Permutation, uncode, Plactic
from schubmult.rings.rc_graph import RCGraph
from schubmult.rings.rc_graph_ring import RCGraphRing


def eg_PQ_from_rc(rc):
    """
    Try to extract (P,Q) = edelman_greene(rc). Many implementations provide
    rc.edelman_greene() returning (P,Q) or similar. Return (P,Q) or raise.
    """
    assert isinstance(rc, RCGraph)
    P, Q = rc.edelman_greene()
    return Plactic(P), Plactic(Q)


def shape_of_tableau(Q):
    """Return the shape (partition) of a tableau-like object Q.
    Accepts objects with `.shape` or sequences of rows (list/tuple)."""
    if Q is None:
        return ()
    if hasattr(Q, "shape"):
        return tuple(Q.shape)
    rows = [len(r) for r in Q]
    while rows and rows[-1] == 0:
        rows.pop()
    return tuple(rows)


def main(argv):
    if len(argv) < 2:
        print("Usage: python check_q_tableau.py N")
        return 2
    try:
        n = int(argv[1])
    except ValueError:
        print("N must be an integer")
        return 2

    length = n - 1
    perms = Permutation.all_permutations(n)
    total_pairs = 0
    failures = []
    summary = Counter()
    rc_ring = RCGraphRing()

    # Exceptions inside the loop are not caught here and will propagate.
    for pu, pv in itertools.product(perms, perms):
        total_pairs += 1
        for a_rc in RCGraph.all_rc_graphs(pu):
            for len1 in range(len(pu.trimcode), n):
                left_rc = a_rc
                if len(a_rc) < len1:
                    left_rc = a_rc.extend(len1 - len(a_rc))
                for b_rc in RCGraph.all_rc_graphs(pv):
                    for len2 in range(len(pv.trimcode), n):
                        right_rc = b_rc
                        if len(b_rc) < len2:
                            right_rc = b_rc.extend(len2 - len(b_rc))

                        # extract EG P/Q for left and shifted-right
                        P_left, Q_left = eg_PQ_from_rc(left_rc)
                        P_right, Q_right = eg_PQ_from_rc(RCGraph(right_rc.shiftup(len(left_rc))))

                        # expected Q is Plactic product of left.Q and shifted right.Q
                        expected_Q = Q_left * Q_right

                        # compute product in RC-graphs using RCGraphRing
                        prod_dict = rc_ring(left_rc) * rc_ring(right_rc)

                        # check every term
                        for rc_res, coeff in prod_dict.items():
                            P_res, Q_res = eg_PQ_from_rc(rc_res)
                            ok_shape = Q_res == expected_Q
                            if not ok_shape:
                                failures.append((pu, pv, rc_res, Q_left, Q_right, Q_res, expected_Q))
                            else:
                                summary["ok"] += 1

    # Report
    print("check_q_tableau report")
    print("----------------------")
    print(f"N = {n}; total permutation pairs considered = {total_pairs}")
    print(f"Total successful matches: {summary['ok']}")
    print(f"Total mismatches found: {len(failures)}")
    if failures:
        print("\nExamples of mismatches (showing up to 10):")
        for ex in failures[:10]:
            pu, pv, rc_res, Q_left, Q_right, got_Q, expected = ex
            print(f"- pair left={pu}, right={pv}")
            print(f"  left.Q = {Q_left}, right.Q = {Q_right}")
            print(f"  result.Q = {got_Q}, expected (Plactic product) = {expected}")
            print(f"  resulting RCGraph repr: {repr(rc_res)}\n")

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))