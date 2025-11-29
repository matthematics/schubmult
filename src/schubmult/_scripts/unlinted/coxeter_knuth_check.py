# check if CoxeterKnuth tensor free algebra is isomoprhic to RCGraphRing
import itertools

from schubmult.ck_ring import CoxeterKnuthRing
from schubmult import PlacticAlgebra
from schubmult import RCGraph
from schubmult import RCGraphRing

from schubmult import FA, ASx, Permutation


def hom(ck_elem):
    """
    Map a CoxeterKnuthRing element to an RCGraphRing element.
    """
    bob = ASx.zero
    for (ck, l), coeff in ck_elem.items():
        bob += coeff * ASx(~ck.perm, l)
    return bob

def hom_rc(rc_elem, shift=0):
    """
    Map an RCGraphRing element to a CoxeterKnuthRing element.
    """
    ck_ring = CoxeterKnuthRing()
    result = ck_ring.zero
    for rc, coeff in rc_elem.items():
        result += coeff * ck_ring((rc.p_tableau, rc.weight_tableau, len(rc)))
    return result

if __name__ == "__main__":
    import sys

    from sympy import pretty_print

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(int(sys.argv[1]))
    rc_ring = RCGraphRing()
    ck_ring = CoxeterKnuthRing()
    tring = ck_ring @ PlacticAlgebra()

    for perm1, perm2 in itertools.product(perms, perms):
        if perm1.inv == 0 or perm2.inv == 0:
            continue
        for len1 in range(len(perm1.trimcode), n):
            for len2 in range(len(perm2.trimcode), n):
                for rc1 in RCGraph.all_rc_graphs(perm1, len1):
                    for rc2 in RCGraph.all_rc_graphs(perm2, len2):
                        # ck1_elem = ck_ring((rc1.p_tableau, len(rc1)))
                        # ck2_elem = ck_ring((rc2.p_tableau, len(rc2)))
                        # mul_elem_ck = ck1_elem * ck2_elem
                        mul_elem_rc = rc_ring(rc1) * rc_ring(rc2)
                        mul_elem_as = hom_rc(rc_ring(rc1)) * hom_rc(rc_ring(rc2))
                        try:
                            assert all(v == 0 for v in (mul_elem_as - hom_rc(mul_elem_rc)).values()), "Mismatch in product for mul_elem"
                        except AssertionError as e:
                            print(e)
                            pretty_print(mul_elem_as)
                            pretty_print(mul_elem_rc)
                            pretty_print(hom_rc(mul_elem_rc))
                            pretty_print(hom_rc(rc_ring(rc1)))
                            pretty_print(hom_rc(rc_ring(rc2)))
                            input()
                        # for rc, coeff in mul_elem_rc.items():
                        #     try:
                        #         assert mul_elem_ck.get((rc.p_tableau, len(rc)), 0) == coeff
                        #     except AssertionError:
                        #         print(f"Mismatch found: {mul_elem_ck.get((rc.p_tableau, len(rc)), 0)} != {coeff}")
                                #print(f"{mul_elem_ax=} {mul_elem_ck=}")
                #         rc_prod = rc_ring(rc1) * rc_ring(rc2)
                #         ck1 = rc1.p_tableau
                #         ck2 = rc2.p_tableau
                #         test1_elem = test_ring(((ck1,len(rc1)), rc1.length_vector))
                #         test2_elem = test_ring(((ck2, len(rc2)), rc2.length_vector))
                #         ck_prod = test1_elem * test2_elem
                #         pretty_print(rc_prod)
                #         pretty_print(ck_prod)

                #         for rc0, coeff in rc_prod.items():
                #             ck0 = rc0.p_tableau
                #             try:
                #                 assert ck_prod.get(((ck0, len(rc0)), rc0.length_vector), 0) == coeff, f"Mismatch in product for {test1_elem} * {test2_elem} at {rc0}"
                #             except AssertionError as e:
                #                 print(e)
                #                 pretty_print(rc_prod)
                #                 pretty_print(ck_prod)
                #                 print(f"{coeff=}")
                #                 print(f"{ck_prod.get(((ck0, len(rc0)), rc0.length_vector), 0)=}")
                #                 pretty_print(rc0)
                #                 pretty_print(rc0.p_tableau)
                #                 input()
