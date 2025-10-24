# check if CoxeterKnuth tensor free algebra is isomoprhic to RCGraphRing
from schubmult.rings.ck_ring import CoxeterKnuthRing
from schubmult.rings.rc_graph_ring import RCGraphRing
from schubmult.rings.rc_graph import RCGraph
from schubmult import Permutation, FA
import itertools

if __name__ == "__main__":
    import sys
    from sympy import pretty_print

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(int(sys.argv[1]))
    ck_ring = CoxeterKnuthRing()
    rc_ring = RCGraphRing()
    test_ring = ck_ring @ FA
    for perm1, perm2 in itertools.product(perms, perms):
        for len1 in range(len(perm1.trimcode), n):
            for len2 in range(len(perm2.trimcode), n):
                for rc1 in RCGraph.all_rc_graphs(perm1, len1):
                    for rc2 in RCGraph.all_rc_graphs(perm2, len2):
                        rc_prod = rc_ring(rc1) * rc_ring(rc2)
                        ck1 = rc1.p_tableau
                        ck2 = rc2.p_tableau
                        test1_elem = test_ring(((ck1,len(rc1)), rc1.length_vector))
                        test2_elem = test_ring(((ck2, len(rc2)), rc2.length_vector))
                        ck_prod = test1_elem * test2_elem
                        pretty_print(rc_prod)
                        pretty_print(ck_prod)

                        for rc0, coeff in rc_prod.items():
                            ck0 = rc0.p_tableau
                            assert ck_prod.get(((ck0, len(rc0)), rc0.length_vector), 0) == coeff, f"Mismatch in product for {test1_elem} * {test2_elem} at {rc0}"
                