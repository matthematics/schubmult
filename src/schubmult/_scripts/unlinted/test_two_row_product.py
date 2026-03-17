from schubmult import *
from schubmult.utils.perm_utils import weak_compositions


if __name__ == "__main__":
    import sys
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.rings.polynomial_algebra import KeyPolyBasis, PolynomialAlgebra, SchubertPolyBasis, Schub
    import itertools

    num_parts = int(sys.argv[1])
    max_part_size = int(sys.argv[2])
    do_assoc = False
    comps = weak_compositions(num_parts, max_part_size)

    r = DualRCGraphRing()

    for comp1, comp2 in itertools.product(comps, repeat=2):
        perm1 = uncode(comp1)
        perm2 = uncode(comp2)
        #perm3 = uncode(comp3)
        # if sum(comp1) == 0 or sum(comp2) == 0 or sum(comp3) == 0:
        #     continue
        #if perm1.is_dominant or perm2.is_dominant or perm3.is_dominant:
        #if perm2.is_dominant or perm3.is_dominant:
        #if perm3.is_dominant:
        #if perm1.is_dominant or perm3.is_dominant:
        # if perm1.is_dominant or perm2.is_dominant:
        #     continue
        
        produ = Schub(perm1, num_parts) * Schub(perm2, num_parts)
        #result = r()
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, num_parts), RCGraph.all_rc_graphs(perm2, num_parts)):
            print("I iz da testing")
            print(rc1)
            print(rc2)
            elem1 = r(rc1)
            elem2 = r(rc2)
            ring_result = elem1 * elem2
            result = next(iter(ring_result))
            assert produ.get((result.perm, num_parts), 0) != 0, f"Failure for comp {comp1}, {comp2}, got {result.perm=}, {produ=}, {result=}, {rc1=} {rc2=}\n{ring_result=}\n{elem1=}\n{elem2=}"
            print("Bimodmod")
            #rc2_base, rc2_grass = rc2.squash_decomp()
            # rc1_base, rc1_grass = rc1.squash_decomp()
            # rc2_base, rc2_grass = rc2.squash_decomp()
            # if rc1_base.perm.inv == 0 or rc2_base.perm.inv == 0:
            #     continue
            # result2 = next(iter(r(rc1_base) *r(rc2.left_squash(rc1_grass))))
            # result2 = r.key_to_rc_graph(result2)
            # assert result == result2, f"Failure for comp {comp1}, {comp2}, got {result=}, {result2=}, {rc1=} {rc2=}"
        #produ = Schub(perm1,2) * Schub(perm2,2) * Schub(perm3,2)
        
        
        #ASSOC
        if do_assoc:
            for comp3 in comps:
                perm3 = uncode(comp3)
                for rc3 in RCGraph.all_rc_graphs(perm3, num_parts):
                    # rc3_base, rc3_grass = rc3.squash_decomp()
                    # if rc3_base.perm.inv == 0:
                    #     continue
                    produ3 = produ * Schub(perm3, num_parts)
                    elem3 = r.from_rc_graph(rc3)
                    # move rc1 across
                    # rc1_base, rc1_grass = rc1.squash_decomp()
                    # # righto_mid = (r(rc1_grass) * r(rc2))
                    # # right_mid_base, right_mid_grass = next(iter(righto_mid)).squash_decomp()
                    # # throughit = r(right_mid_grass) * r(rc3)
                    # # throughit_base, throughit_grass = next(iter(throughit)).squash_decomp()
                    # # righto = (r(right_mid_base) * r(throughit_base)) * r(throughit_grass)
                    # # righto_o_base, righto_o_grass = next(iter(righto)).squash_decomp()
                    # # result2 = next(iter((r(rc1_base) * r(righto_o_base))* r(righto_o_grass)))
                    # righto2 = r(rc1_grass) * (r(rc2) * r(rc3))
                    # #assert righto == righto2, f"Failure for comp {comp1}, {comp2}, {comp3}, got {righto=}, {righto2=}"
                    # righto_base, righto_grass = next(iter(righto2)).squash_decomp()
                    # result1 = next(iter((r(rc1_base)*r(righto_base)) * r(righto_grass)))
                    result = elem1 * (elem2 * elem3)
                    result1 = next(iter(result))
                    result1 = r.key_to_rc_graph(result1)
                    assert produ3.get((result1.perm, num_parts), 0) != 0, f"Failure for comp {comp1}, {comp2}, {comp3}, got {result1.perm=}, {produ3=}, {result1=}\n{(r(rc1) * r(rc2))} * {r(rc3)}\n{r(rc1)} * {(r(rc2) * r(rc3))}"
                    result2 = (elem1 * elem2) * elem3
                    result2 = next(iter(result2))
                    result2 = r.key_to_rc_graph(result2)
                    assert produ3.get((result2.perm, num_parts), 0) != 0, f"Failure for comp {comp1}, {comp2}, {comp3}, got {result2.perm=}, {produ3=}, {result2=}\n{(r(rc1) * r(rc2))} * {r(rc3)}\n{r(rc1)} * {(r(rc2) * r(rc3))}"
                    #result2 = next(iter(r(rc1_base) * (r(rc1_grass) * (r(rc2) * r(rc3)))))
                    
                    #assert produ3.get((result1.perm, num_parts), 0) != 0, f"Failure for comp {comp1}, {comp2}, {comp3}, got {result1.perm=}, {produ3=}, {result1=}\n{(r(rc1) * r(rc2))} * {r(rc3)}\n{r(rc1)} * {(r(rc2) * r(rc3))}"
                    
                    # left_factor = next(iter(r(rc1) * r(rc2)))
                    # left_factor_base, left_factor_grass = left_factor.squash_decomp()
                    # right_factor = rc3.left_squash(left_factor_grass)
                    # right_factor_base, right_factor_grass = right_factor.squash_decomp()
                    # result2 = next(iter((r(left_factor_base) * r(right_factor_base)) * r(right_factor_grass)))
                    # assert produ3.get((result2.perm, num_parts), 0) != 0, f"Failure for comp {comp1}, {comp2}, {comp3}, got {result2.perm=}, {produ3=}, {result2=}\n{(r(rc1) * r(rc2))} * {r(rc3)}\n{r(rc1)} * {(r(rc2) * r(rc3))}"
                    # assert result1 == result2, f"Failure for comp {comp1}, {comp2}, {comp3}, got {result1=}, {result2=}, {produ3=}\n{(r(rc1) * r(rc2))} * {r(rc3)}\n{r(rc1)} * {(r(rc2) * r(rc3))}"
                    print("Success for comp", comp1, comp2, comp3)
        print("Full success for comp", comp1, comp2)
        #         bank1 = (r(rc1_grass)*r(rc2))
        #         bank1_base, bank1_grass = next(iter(bank1)).squash_decomp()
        #         bank1_right = r(bank1_grass) * r(rc3)
        #         bank1_right_base,  bank1_right_grass = next(iter(bank1_right)).squash_decomp()
        #         # bank1_right_base, bank1_right_grass = next(iter(bank1_right)).squash_decomp()
        #         result2 = next(iter(((r(rc1_base) * r(bank1_base))*r(bank1_right_base))*r(bank1_right_grass)))
        #         assert produ3.get((result2.perm, num_parts), 0) != 0, f"Failure for comp {comp1}, {comp2}, {comp3}, got {result2.perm=}, {produ3=}, {result2=}\n{(r(rc1) * r(rc2))} * {r(rc3)}\n{r(rc1)} * {(r(rc2) * r(rc3))}"
        #         assert result1 == result2, f"Failure for comp {comp1}, {comp2}, {comp3}, got {result1=}, {result2=}\n{(r(rc1) * r(rc2))} * {r(rc3)}\n{r(rc1)} * {(r(rc2) * r(rc3))}"
        #         #assert produ3.get((result2.perm, num_parts), 0) != 0, f"Failure for comp {comp1}, {comp2}, {comp3}, got {result2.perm=}, {produ3=}, {result2=}\n{(r(rc1) * r(rc2))} * {r(rc3)}\n{r(rc1)} * {(r(rc2) * r(rc3))}"
        #         #assert result1 == result2, f"Failure for comp {comp1}, {comp2}, {comp3}, got {result1=}, {result2=}\n{(r(rc1) * r(rc2))} * {r(rc3)}\n{r(rc1)} * {(r(rc2) * r(rc3))}"
        # # for rc1, rc2, rc3 in itertools.product(RCGraph.all_rc_graphs(perm1, num_parts), RCGraph.all_rc_graphs(perm2, num_parts), RCGraph.all_rc_graphs(perm3, num_parts)):
        # #     rc1_base, rc1_grass = rc1.squash_decomp()
        # #     middle = rc2.left_squash(rc1_grass)
        # #     middle_base, middle_grass = middle.squash_decomp()
        # #     right = rc3.left_squash(middle_grass)
        # #     right_base, right_grass = right.squash_decomp()

        # #     result = next(iter(((r(rc1_base) * r(middle_base)) * r(right_base)) * r(right_grass)))
        # #     result2 = next(iter(((r(rc1_base) * (r(middle_base) * r(right_base))) * r(right_grass))))
        # #     #assert produ.get((result.perm, 2), 0) != 0, f"Failure for comp {comp1}, {comp2}, got {result.perm=}, {produ=}, {result=}"
        # #     #assert result2.almosteq(result), f"Failure for comp {comp1}, {comp2}, got {result.perm=}, {produ=}, {result=}, {result2=}"
        # #     assert result2 == result, f"Failure for comp {comp1}, {comp2}, {comp3}, got {result.perm=}, {result=}, {result2=}\n{(r(rc1) * r(rc2))} * {r(rc3)}\n{r(rc1)} * {(r(rc2) * r(rc3))}"
        # #     assert produ.get((result.perm, 2), 0) != 0, f"Failure for comp {comp1}, {comp2}, {comp3}, got {result.perm=}, {produ=}, {result=}\n{(r(rc1) * r(rc2))} * {r(rc3)}\n{r(rc1)} * {(r(rc2) * r(rc3))}"
        # #     print("Success for comp", comp1, comp2, comp3)