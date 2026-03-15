from schubmult import *
from schubmult.utils.perm_utils import weak_compositions


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
        produ = Schub(perm1,2) * Schub(perm2,2)
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, num_parts), RCGraph.all_rc_graphs(perm2, num_parts)):
            print("I iz da testing")
            print(rc1)
            print(rc2)
            two_rc, three_rc = rc1.squash_decomp()
            two_rc2, three_rc2 = rc2.squash_decomp()
            #bimodule tensor product
            result = next(iter((r(two_rc.resize(2))*r(two_rc2.resize(2))).resize(3)*(r(three_rc)*r(three_rc2))))
            #next(iter(r(two_rc.resize(2)) * r(two_rc2.resize(2)))).resize(3).squash_product(three_rc).squash_product(three_rc2)
            middle = (three_rc.left_squash(two_rc2)).squash_product(three_rc2)
            middle_base, middle_grass = middle.squash_decomp()
            # tensor_part = r(two_rc.resize(2)) * r(middle_base.resize(2))
            # tensor_result = next(iter(tensor_part)).resize(3)
            # tensor_result_hw, raise_seq = tensor_result.to_highest_weight()
            # tensor_hw_base, tensor_hw_grass = tensor_result_hw.squash_decomp()
            # right_grass_part = tensor_hw_grass.squash_product(middle_grass)
            # new_tensor = CrystalGraphTensor(tensor_hw_base, right_grass_part).reverse_raise_seq(raise_seq)
            #result = tensor_base.squash_product(right_grass_part)
            #result = new_tensor.factors[0].squash_product(new_tensor.factors[1])
            #hw.factors[0].squash_product(hw.factors[1]).reverse_raise_seq(raise_seq)
            #result_base, result_grass = next(iter((r(two_rc.resize(2)) * r(middle_base.resize(2))))).resize(3).squash_decomp()
            
            result = next(iter((r(two_rc.resize(2)) * r(middle_base.resize(2))))).resize(3).squash_product(middle_grass)
            # # result = next(iter(r(rc1) * r(rc2)))
            assert produ.get((result.perm, 2), 0) != 0, f"Failure for comp {comp1}, {comp2}, got {result.perm=}, {produ=}, {result=}, {rc1=} {rc2=}"
            print("Success for comp", comp1, comp2)