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

def grass_perm(shape, length):
    sputnik = list(reversed(shape))
    while len(sputnik) < length:
        sputnik = [0, *sputnik]
    return uncode(tuple(sputnik))
from schubmult import Sx
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
    from sympy import pretty_print
    for pu, pv in itertools.combinations(perms, 2):
        for left_rc in RCGraph.all_rc_graphs(pu):
            for right_rc in RCGraph.all_rc_graphs(pv):
                #left_rc_hw, _ = left_rc.to_highest_weight()
                P_left, Q_left = eg_PQ_from_rc(left_rc)
                P_right, Q_right = eg_PQ_from_rc(right_rc)

                lefto_set = left_rc.prod_with_rc(RCGraph([()] * (len(right_rc.perm))))
                for lefto in lefto_set:
                    while len(lefto) > len(left_rc) + len(right_rc):
                        lefto = lefto.zero_out_last_row()
                    P_lefto, Q_lefto = eg_PQ_from_rc(lefto)
                    assert P_lefto.shape == P_left.shape
                    # then it's right times left
                # def calc_the_baby(left, right):
                #     len0 = max(len(left.perm.trimcode), len(right.perm.trimcode)) 
                #     left_rc_hw = left
                #     right_rc_hw = right
                #     lrc = (~(left_rc_hw)).normalize()
                #     rrc = (~(right_rc_hw)).normalize()
                #     tprods =   rc_ring(lrc)* rc_ring(rrc)
                #     tprodst = rc_ring.zero
                #     for rc1, coeff1 in tprods.items():
                #         pain = (~rc1)
                        
                #         if len(pain) < max(len0, len(pain.perm.trimcode)):

                #             pain = pain.resize(max(len0, len(pain.perm.trimcode))+5)
                #             while len(pain) > len0 and len(pain[-1]) == 0:
                #                 pain = pain.zero_out_last_row()
                #         if len(pain) == len0:# and pain.length_vector == tuple([a+b for a,b in zip_longest(left_rc_hw.length_vector, right_rc_hw.length_vector, fillvalue=0)]):
                #             tprodst += coeff1 * rc_ring(pain)
                #     return tprodst
                
                # pretty_print(left_rc)
                # print("tensor")
                # pretty_print(right_rc)
                # print("=====> product")
                # tprodst = calc_the_baby(left_rc, right_rc)
                # pretty_print(tprodst)
                # tprodst2 = calc_the_baby(right_rc, left_rc)
                # pretty_print(tprodst2)
                # print("===================")

    #                     # extract EG P/Q for left and right (may be plain table objects)
    #                     P_left, Q_left = eg_PQ_from_rc(left_rc)
    #                     P_right, Q_right = eg_PQ_from_rc(right_rc)

    #                     assert P_left.perm == ~left_rc.perm
    #                     assert P_right.perm == ~right_rc.perm

                        
    #                     # normalize Q's to Plactic instances (accept tuple-of-rows or Plactic

    #                     # shift entries of right-Q up by len1 (add len1 to every entry)

    #                     # compute expected Q in the "standard" domain if semistandard ordering is reversed
    #                     # expected_Q = Plactic.from_standard(
    #                     #     Q_left.to_standard(reverse_semistandard=True)
    #                     #            * Q_right.to_standard(reverse_semistandard=True),
    #                     #     reverse_semistandard=True,
    #                     # )
    #                     # expected_Q = Q_left * Q_right

    #                     # # expected Q is Plactic product of left.Q and shifted right.Q
    #                     # #expected_Q = Q_left * Q_right
    #                     # print(f"{left_rc=}")
    #                     # print(f"{P_left=}")
    #                     # print(f"{Q_left=}")
    #                     # print(f"{right_rc=}")
    #                     # print(f"{P_right=}")
    #                     # print(f"{Q_right=}")
    #                     # # compute product in RC-graphs using RCGraphRing
    #                     prod_dict = rc_ring(left_rc) * rc_ring(right_rc)
    #                     # print(f"{left_rc.perm.trimcode=}")
    #                     # print(f"{right_rc.perm.trimcode=}")
    #                     from schubmult import Sx
    #                     # check every term
    #                     for rc_res, coeff in prod_dict.items():
    #                         if coeff == 0:
    #                             continue
    #                         # if len(rc_res.perm.descents()) > 1:
    #                         #     continue  # skip non-Grassmannian results for now
    #                         P_res, Q_res = eg_PQ_from_rc(rc_res)
    #                         assert Sx(grass_perm(P_res.shape,len1+len2)).coproduct(*list(range(1, len1+1))).get((grass_perm(P_left.shape,len1),grass_perm(P_right.shape,len2)),0) != 0#, f"Product P-shape mismatch {dict(Sx(grass_perm(P_res.shape)).coproduct(*list(range(len1+1,len1+len2+1))))} {grass_perm(P_res.shape)=} {grass_perm(P_left.shape)=} {grass_perm(P_right.shape)=} {rc_res.perm.trimcode=} {left_rc.perm.trimcode=} {right_rc.perm.trimcode=}"
    #                     print("Success")
    #                         # ok_shape = (Q_res == expected_Q)
    #                         # if not ok_shape:
    #                         #     failures.append((pu, pv, rc_res, Q_left, Q_right, Q_res, expected_Q))
    #                         #     print("Fail")
    #                         #     print(f"{rc_res.perm.trimcode=}")
    #                         #     print(f"{P_res=}")
    #                         #     print(f"{Q_res=}")
    #                         #     print(f"{expected_Q=}")
    #                         # else:
    #                         #     summary["ok"] += 1
    #                         #     print("Ok")

    # # Report
    # # print("check_q_tableau report")
    # # print("----------------------")
    # # print(f"N = {n}; total permutation pairs considered = {total_pairs}")
    # # print(f"Total successful matches: {summary['ok']}")
    # # print(f"Total mismatches found: {len(failures)}")
    # # if failures:
    # #     print("\nExamples of mismatches (showing up to 10):")
    # #     for ex in failures[:10]:
    # #         pu, pv, rc_res, Q_left, Q_right, got_Q, expected = ex
    # #         print(f"- pair left={pu}, right={pv}")
    # #         print(f"  left.Q = {Q_left}, right.Q = {Q_right}")
    # #         print(f"  result.Q = {got_Q}, expected (Plactic product) = {expected}")
    # #         print(f"  resulting RCGraph repr: {repr(rc_res)}\n")

    # # return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))