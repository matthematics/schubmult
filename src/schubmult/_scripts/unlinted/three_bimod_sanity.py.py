from schubmult import *
from schubmult.utils.perm_utils import weak_compositions
from schubmult.utils.tuple_utils import pad_tuple

def is_full_grass(perm, n):
    return perm.inv == 0 or perm.descents() == {n - 1}


if __name__ == "__main__":
    import sys
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.rings.polynomial_algebra import KeyPolyBasis, PolynomialAlgebra, SchubertPolyBasis, Schub
    import itertools

    num_parts = 3
    max_part_size = int(sys.argv[1])

    comps = weak_compositions(num_parts, max_part_size)

    r = DualRCGraphRing()

    
    for comp1, comp2, comp3 in itertools.product(comps, repeat=3):
        perm1 = uncode(comp1)
        perm2 = uncode(comp2)
        perm3 = uncode(comp3)
        if len(perm2) > num_parts:
            continue
        if not is_full_grass(perm1, num_parts) or not is_full_grass(perm3, num_parts) or is_full_grass(perm2, num_parts):
            continue
        
        for rc1, rc2, rc3 in itertools.product(RCGraph.all_rc_graphs(perm1, num_parts), RCGraph.all_hw_rcs(perm2, num_parts), RCGraph.all_rc_graphs(perm3, num_parts)):
            result1 = rc2.left_squash(rc1).squash_product(rc3)
            result2 = rc2.left_squash(rc1.squash_product(rc3))
            assert result1 == result2, f"Failure for comp {comp1}, {comp2}, {comp3}, {result1=}, {result2=}"

        print("Success for comp", comp1, comp2, comp3)