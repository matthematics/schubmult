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
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, num_parts), RCGraph.all_rc_graphs(perm2, num_parts)):
            print("I iz da testing")
            print(rc1)
            print(rc2)
            rc1_base, rc1_grass = rc1.squash_decomp()
            if rc1_base.perm.inv == 0:
                continue
            rc2_base, rc2_grass = rc2.squash_decomp()
            if rc2_base.perm.inv == 0:
                continue
            if rc1_grass.perm.inv != 0:
                right_factor = rc2.left_squash(rc1_grass)
            else:
                right_factor = rc2
            middle_base, rightmost_grass = right_factor.squash_decomp()
            if middle_base.perm.inv == 0:
                result = rc1_base.squash_product(right_factor)
            elif rightmost_grass.perm.inv == 0:
                result = next(iter(r(rc1_base.resize(2)) * r(middle_base.resize(2)))).resize(3)
            else:
                two_result = next(iter((r(rc1_base.resize(2)) * r(middle_base.resize(2))).resize(3)))
                result = two_result.squash_product(rightmost_grass)
            # rc2_base, rc2_grass = rc2.squash_decomp()
            # print(f"{rc1_base=}, {rc1_grass=}")
            # print(f"{rc2_base=}, {rc2_grass=}")
            # if rc2_base.perm.inv == 0:
            #     result = rc1.squash_product(rc2_grass)
            # else:
            #     right_fact = rc1_grass.left_squash(rc2)
            #     middle_base, rightmost_grass = right_fact.squash_decomp()
            #     #rightmost_grass = middle_grass.squash_product(rc2_grass)
            #     # rc1_base, middle_base, middle_grass, rc2_grass
            #     assert len(rc1_base.perm) <= 3
            #     assert len(middle_base.perm) <= 3
                
            #     print("Yay hoory")
            #     # continue

            #     #rc1_base and right_fact
            #     # tensor crystal
            #     u_rc, v_rc = rc1_base, middle_base
            #     two_result = None
            #     if u_rc.perm.is_dominant and v_rc.perm.is_dominant:
            #         two_result = RCGraph.principal_rc(uncode([a + b for a, b in zip(pad_tuple(u_rc.perm.trimcode, num_parts), pad_tuple(v_rc.perm.trimcode, num_parts))]), 3)
            #     elif not v_rc.perm.is_dominant:
            #         tensor = CrystalGraphTensor(u_rc, v_rc)
            #         tensor_hw, raise_seq = tensor.to_highest_weight()
            #         if tensor == tensor_hw and tensor.is_lowest_weight:
            #             two_result = RCGraph.principal_rc(uncode(tensor_hw.crystal_weight), 3)
            #         else:
            #             two_result_hw = next(iter(RCGraph.all_hw_rcs(uncode(tuple(reversed(tensor_hw.crystal_weight))), 3)))
            #             two_result = two_result_hw.reverse_raise_seq(raise_seq).to_rc_graph()
            #     else:
            #         tensor = CrystalGraphTensor(AntiRCGraph.from_rc_graph(u_rc), AntiRCGraph.from_rc_graph(v_rc))
            #         tensor_hw, raise_seq = tensor.to_highest_weight()
            #         if tensor == tensor_hw:
            #             two_result = RCGraph.principal_rc(uncode(tuple(reversed(tensor_hw.crystal_weight))), 3)
            #         else:
            #             two_result_hw = AntiRCGraph.from_rc_graph(next(iter(RCGraph.all_hw_rcs(uncode(tuple((tensor_hw.crystal_weight))), 3))))
            #             two_result = two_result_hw.reverse_raise_seq(raise_seq).to_rc_graph()

            #     print(f"{u_rc=}")
            #     print(f"{v_rc=}")
            #     print(f"{two_result=}")
            #     print(f"{rightmost_grass=}")
            #     #two_result_base, two_result_grass = two_result.squash_decomp()
            #     #result = two_result.squash_product(rightmost_grass)
            #     tensor2 = CrystalGraphTensor(two_result, rightmost_grass)
            #     tensor2_hw, raise_seq2 = tensor2.to_highest_weight()
            #     # if tensor2 == tensor2_hw and tensor2.is_lowest_weight:
            #     #     result = RCGraph.principal_rc(uncode(tensor2_hw.crystal_weight), 3)
            #     # else:
            #     #     result_hw = next(iter(RCGraph.all_hw_rcs(uncode(tuple(reversed(tensor2_hw.crystal_weight))), 3)))
            #     #     result = result_hw.reverse_raise_seq(raise_seq2)
            #     result = tensor2_hw.factors[0].squash_product(tensor2_hw.factors[1]).reverse_raise_seq(raise_seq2)
            # two_rc2, three_rc2 = rc2.squash_decomp()
            # #bimodule tensor product
            # # result = next(iter((r(two_rc.resize(2))*r(two_rc2.resize(2))).resize(3)*(r(three_rc)*r(three_rc2))))
            # #next(iter(r(two_rc.resize(2)) * r(two_rc2.resize(2)))).resize(3).squash_product(three_rc).squash_product(three_rc2)
            # middle = (three_rc.left_squash(two_rc2))#.squash_product(three_rc2)
            # middle_base, middle_grass = middle.squash_decomp()
            # right_grass = three_rc2.left_squash(middle_grass)
            #middle_grass.lesquash(three_rc2)
            # two_rc2 = two_rc.resize(2)
            # middle_base2 = middle_base.resize(2)
            # # assert len(two_rc2.perm) <= 3
            # # assert len(middle_base2.perm) <= 3
            # two_rc2_base, two_rc2_grass = two_rc2.squash_decomp()
            # middle_base2_middle = two_rc2_grass.left_squash(middle_base2)
            # middle_base2_middle_base, middle_base2_middle_grass = middle_base2_middle.squash_decomp()
            # left_factor2 = next(iter(r(two_rc2_base) * r(middle_base2_middle_base))).squash_product(middle_base2_middle_grass)
            # left_factor = left_factor2.resize(3)
            # result = left_factor.squash_product(middle_grass)
            # middle_base2_base, middle_base2_grass = middle_base2.squash_decomp()
            # middle_base2_middle = two_rc2_base.left_squash(middle_base2_base)
            # middle_base2_middle_base, middle_base2_middle_grass = middle_base2_middle.squash_decomp()
            # left_factor2 = next(iter(r(two_rc2) * r(middle_base2)))
            # left_factor = left_factor2.resize(3)
            # left_base, left_grass = left_factor.squash_decomp()
            # rightmost_grass = left_grass.left_squash(middle_grass)
            # # two_rc_base, two_rc_grass = two_rc2.squash_decomp()
            # result = left_base.squash_product(rightmost_grass)
            # middle_base2 = middle_base.resize(2)
            # middle2 = two_rc_grass.left_squash(middle_base2)
            # right_base2, right_grass2 = middle2.squash_decomp()

            #righ

            #result  = next(iter((r(two_rc.resize(2)) * r(middle_base.resize(2))).resize(3))).squash_product(middle_grass)
            # tensor_result = next(iter(tensor_part)).resize(3)
            # tensor_result_hw, raise_seq = tensor_result.to_highest_weight()
            # tensor_hw_base, tensor_hw_grass = tensor_result_hw.squash_decomp()
            # right_grass_part = tensor_hw_grass.squash_product(middle_grass)
            # new_tensor = CrystalGraphTensor(tensor_hw_base, right_grass_part).reverse_raise_seq(raise_seq)
            #result = tensor_base.squash_product(right_grass_part)
            #result = new_tensor.factors[0].squash_product(new_tensor.factors[1])
            #hw.factors[0].squash_product(hw.factors[1]).reverse_raise_seq(raise_seq)
            #result_base, result_grass = next(iter((r(two_rc.resize(2)) * r(middle_base.resize(2))))).resize(3).squash_decomp()
            
            #result = next(iter((r(two_rc.resize(2)) * r(middle_base.resize(2))))).resize(3).squash_product(right_grass)
            # # result = next(iter(r(rc1) * r(rc2)))
            assert produ.get((result.perm, 3), 0) != 0, f"Failure for comp {comp1}, {comp2}, got {result.perm=}, {produ=}, {result=}, {rc1=} {rc2=}"
            print("Success for comp", comp1, comp2)