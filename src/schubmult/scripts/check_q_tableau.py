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

from schubmult.perm_lib import Permutation, uncode
from schubmult.rings.rc_graph import RCGraph
from schubmult.rings.rc_graph_ring import RCGraphRing
from schubmult.rings.plactic import Plactic


def eg_PQ_from_rc(rc):
    """
    Try to extract (P,Q) = edelman_greene(rc). Many implementations provide
    rc.edelman_greene() returning (P,Q) or similar. Return (P,Q) or raise.
    """
    assert isinstance(rc, RCGraph)
    return rc.edelman_greene()


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
        if len(pu.descents()) > 1 or len(pv.descents()) > 1:
            continue  # limit to Grassmannian permutations for now
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

                        # extract EG P/Q for left and right (may be plain table objects)
                        P_left, Q_left = eg_PQ_from_rc(left_rc)
                        P_right, Q_right = eg_PQ_from_rc(right_rc)

                        assert P_left.perm == ~left_rc.perm
                        assert P_right.perm == ~right_rc.perm

                        # normalize Q's to Plactic instances (accept tuple-of-rows or Plactic)
                        if not isinstance(Q_left, Plactic):
                            Q_left = Plactic.from_word(Q_left)
                        if not isinstance(Q_right, Plactic):
                            Q_right = Plactic.from_word(Q_right)

                        # shift entries of right-Q up by len1 (add len1 to every entry)
                        Q_right = Plactic(Plactic._remap_word(Q_right._word, lambda x: int(x) + len1))

                        # compute expected Q in the "standard" domain if semistandard ordering is reversed
                        # expected_Q = Plactic.from_standard(
                        #     Q_left.to_standard(reverse_semistandard=True)
                        #            * Q_right.to_standard(reverse_semistandard=True),
                        #     reverse_semistandard=True,
                        # )
                        expected_Q = Q_left * Q_right

                        # expected Q is Plactic product of left.Q and shifted right.Q
                        #expected_Q = Q_left * Q_right
                        print(f"{P_left=}")
                        print(f"{Q_left=}")
                        print(f"{P_right=}")
                        print(f"{Q_right=}")
                        # compute product in RC-graphs using RCGraphRing
                        prod_dict = rc_ring(left_rc) * rc_ring(right_rc)
                        print(f"{left_rc.perm.trimcode=}")
                        print(f"{right_rc.perm.trimcode=}")
                        
                        # check every term
                        for rc_res, coeff in prod_dict.items():
                            if len(rc_res.perm.descents()) > 1:
                                continue  # skip non-Grassmannian results for now
                            P_res, Q_res = eg_PQ_from_rc(rc_res)
                            ok_shape = (Q_res == expected_Q)
                            if not ok_shape:
                                failures.append((pu, pv, rc_res, Q_left, Q_right, Q_res, expected_Q))
                                print("Fail")
                                print(f"{rc_res.perm.trimcode=}")
                                print(f"{P_res=}")
                                print(f"{Q_res=}")
                                print(f"{expected_Q=}")
                            else:
                                summary["ok"] += 1
                                print("Ok")

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