from schubmult import *
from schubmult.utils.perm_utils import weak_compositions
from schubmult.utils.tuple_utils import pad_tuple


if __name__ == "__main__":
    import sys
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.rings.polynomial_algebra import KeyPolyBasis, PolynomialAlgebra, SchubertPolyBasis, Schub
    import itertools

    num_parts = 3
    max_part_size = int(sys.argv[1])

    comps = weak_compositions(num_parts, max_part_size)

    r = DualRCGraphRing()

    
    for comp1, comp2 in itertools.product(comps, repeat=2):
        perm1 = uncode(comp1)
        perm2 = uncode(comp2)
        produ = Schub(perm1,3) * Schub(perm2,3)
        #for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, num_parts), RCGraph.all_rc_graphs(perm2, num_parts)):
        for rc1, rc2 in itertools.product(RCGraph.all_hw_rcs(perm1, num_parts), RCGraph.all_hw_rcs(perm2, num_parts)):
            print(rc1)
            print(rc2)
            
            elem1 = r(rc1)
            elem2 = r(rc2)
            result = next(iter(elem1 * elem2))

            assert produ.get((result.perm, 3), 0) != 0, f"Failure for comp {comp1}, {comp2}, got {result.perm=}, {produ=}, {result=}, {repr(rc1)=} {repr(rc2)=}, {result=}\n{elem1=}\n{elem2=}"
            print("Success for comp", comp1, comp2)