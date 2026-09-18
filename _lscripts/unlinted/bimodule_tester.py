from schubmult import *
from schubmult.symbolic import *
from schubmult.utils.perm_utils import weak_compositions


if __name__ == "__main__":
    import sys
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.rings.polynomial_algebra import KeyPolyBasis, PolynomialAlgebra, SchubertPolyBasis, Schub
    import itertools

    KeyPoly = PolynomialAlgebra(KeyPolyBasis(Sx.genset))

    num_parts = int(sys.argv[1])
    max_part_size = int(sys.argv[2])

    comps = weak_compositions(num_parts, max_part_size)

    r = DualRCGraphRing()

    hw_seen = {}
    grass_hw_seen = set()
    for comp1, comp2 in itertools.product(comps, repeat=2):
        
        perm1 = uncode(comp1)
        perm2 = uncode(comp2)
        if len(perm2) > num_parts:
            continue
        if len(perm1) < num_parts + 1:
            continue
        print(f"Testing comp {comp1}, {comp2}...")
        # if len(perm2) < num_parts + 1 or len(perm2) < num_parts + 1:
        #     continue
        #decomps = {}
        result = 0
        hw_seen = set()
        for rc1 in RCGraph.all_rc_graphs(perm1, num_parts):
             for rc2 in RCGraph.all_rc_graphs(perm2, num_parts):
                base, grass =  rc1.squash_decomp()
                tensor_hw = CrystalGraphTensor(base, grass).to_highest_weight()[0]
                switcheroo = tensor_hw.factors[1].left_squash(rc2)
                base2, grass2 = switcheroo.squash_decomp()
                tensor_hw2 = CrystalGraphTensor(base2, grass2).to_highest_weight()[0]                
                if (base, tensor_hw2) not in hw_seen:
                    hw_seen.add((base, tensor_hw2))
                result += (r(base)@r(tensor_hw2.factors[0])@r(tensor_hw2.factors[1]))@ KeyPoly.from_expr(sum([expand_seq(tens.crystal_weight, Sx.genset) for tens in CrystalGraphTensor(base, tensor_hw2).full_crystal]), length=num_parts)
        assert all(v >= 0 for v in result.values()), f"Failure for comp {comp1}, {comp2}, result {result}"
        # for rc in RCGraph.all_hw_rcs(perm2, num_parts):
        #     base, grass =  rc.squash_decomp()
        #     result += Schub.from_expr(Sx(perm1).expand() * r.from_dict(base @ r(grass)
        #     #base @ KeyPoly.from_expr(Sx(perm1).expand() * grass.polyvalue(Sx.genset), length=num_parts)
        # new_result = 0
        # new_new_result = 0
        # hw_seen = set()
        # for ((perm, _), grass_rc), coeff in result.items():
        #     for new_rc in RCGraph.all_rc_graphs(perm, num_parts):
        #         base, grass = new_rc.squash_decomp()
        #         tensor_hw = CrystalGraphTensor(base, grass.squash_product(grass_rc)).to_highest_weight()[0]
        #         if tensor_hw not in hw_seen:
        #             hw_seen.add(tensor_hw)
        #         new_result += coeff * r(tensor_hw.factors[0]) @ r(tensor_hw.factors[1])
        #         # base_hw = base.to_highest_weight()[0]
        #         # if base_hw not in hw_seen:
        #         #     hw_seen.add(base_hw)
        #         # new_result += coeff * r(base_hw) @ (r(grass) * r(grass_rc))
        #         #new_new_result += coeff * (r(base) * r(grass) * r(grass_rc)
        # # for (rc, grass_perm), coeff in new_result.items():
        # #     new_new_result += coeff * r(rc) @ r.schub(grass_perm)
        # assert all(v >= 0 for v in new_result.values()), f"Failure for comp {comp1}, {comp2}, result {new_result}"
        # # prodo = Sx(perm1) * Sx(perm2)
        # # proff = Sx.from_expr(new_new_result.polyvalue(Sx.genset))
        # # assert prodo == proff, f"Failure for comp {comp1}, {comp2}, new_new_result {new_new_result}, prodo {prodo}, proff {proff}"
        # #assert all(prodo.get(rc.perm, 0) == v for rc, v in new_new_result.items()), f"Failure for comp {comp1}, {comp2}, new_new_result {new_new_result}, prodo {prodo}"
        
        #     # tensor = CrystalGraphTensor(base, grass)
            # tensor_hw = tensor.to_highest_weight()[0]
            # hw = rc.to_highest_weight()[0]
            # if hw not in hw_seen:
            #     hw_seen[hw] = tensor_hw
            # else:
            #     assert tensor_hw == hw_seen[hw], f"Failure for comp {comp}, base {base}, grass {grass}, hw {hw}, tensor_hw {tensor_hw}"
            # grass_hw_seen.add(grass)
                
            #decomps[base] = decomps.get(base, r.zero) + r(grass)
    # results = {}
    # for hw1, hw2 in itertools.product(hw_seen.keys(), grass_hw_seen):
    #     t1 = hw_seen[hw1]
    #     #t2 = hw_seen[hw2]
        # results[(hw1, hw2)] = results.get((hw1, hw2), 0)
        # for tt1 in t1.full_crystal:
        #     for tt2 in t2.full_crystal:
        #         # we have tt1.factors[0] @ tt1.factors[1] * tt2.factors[0] @ tt2.factors[1]
        #         middle = tt1.factors[1].left_squash(tt2.factors[0])
        #         base, grass2 = middle.squash_decomp()
        #         right = grass2.squash_product(tt2.factors[1])
                # #results[(hw1.perm, hw2.perm)] = results.get((hw1.perm, hw2.perm), 0) + KeyPoly.from_expr(tt1.factors[0].polyvalue(Sx.genset) * tt2.factors[0].polyvalue(Sx.genset), length=num_parts) @ (r(tt1.factors[1])*r(tt2.factors[1]))
                
                # #results[(hw1, hw2)] += r(tt1.factors[0]) @ r(base) @ right
                
    #     results[(hw1, hw2)] += r(t1.factors[0]) @ (r(t1.factors[1])*r(hw2))
    #             #Schub.from_expr(tt1.factors[0].polyvalue(Sx.genset) * tt2.factors[0].polyvalue(Sx.genset), length=num_parts) @ (r(tt1.factors[1])*r(tt2.factors[1]))
    #             #+ KeyPoly.from_expr(tt1.factors[0].polyvalue(Sx.genset) * tt1.factors[1].polyvalue(Sx.genset), length=num_parts)@ KeyPoly.from_expr((r(tt2.factors[0])*r(tt2.factors[1])).polyvalue(Sx.genset), length=num_parts)
    # int_results = {}
    # for (hw1, hw2), result in results.items():
    #     for (key1, key2, key3), coeff in result.items():
            
    #         #int_results[(hw1.perm, hw2.perm)] = int_results.get((hw1.perm, hw2.perm), 0) + Schub.from_expr(key1.polyvalue(Sx.genset) * key2.polyvalue(Sx.genset), length=num_parts) @ r(key3)
    #         #int_results[(hw1.perm, hw2.perm)] = int_results.get((hw1.perm, hw2.perm), 0) + coeff * KeyPoly.from_expr(key1.polyvalue(Sx.genset)key3.polyvalue(Sx.genset), length=num_parts)
    #         int_results[(hw1.perm, hw2)] = int_results.get((hw1.perm, hw2), 0) + coeff * Schub.from_expr(key1.polyvalue(Sx.genset)*key2.polyvalue(Sx.genset), length=num_parts) @ r(key3)
    #         #, length=num_parts) @ KeyPoly.from_expr(key3.polyvalue(Sx.genset), length=num_parts)
    # #full_results = {
    # full_results = {}
    # for (perm1, perm2), result in int_results.items():
    #     for ((perm, length), grass_rc), coeff in result.items():
    #         for new_rc in RCGraph.all_rc_graphs(perm, num_parts):
    #             base, grass2_rc = new_rc.squash_decomp()
    #             #fabrika = r(base) * r(grass2_rc) * r(grass_rc)
    #             full_results[(perm1, perm2)] = full_results.get((perm1, perm2), 0) + coeff * KeyPoly.from_expr(base.polyvalue(Sx.genset)) @ (r(grass2_rc) * r(grass_rc))
    #                 #base) @ (r(grass2_rc) * r(grass_rc))
    #     # for (base1, base2, grass2), coeff in result.items():
    #     #     # for new_rc in RCGraph.all_rc_graphs(perm, num_parts):
    #     #     #     base, grass = new_rc.squash_decomp()
    #     #     full_results[(hw1, hw2)] = full_results.get((hw1, hw2), 0) + coeff * Schub.from_expr(base1.polyvalue(Sx.genset) * base2.polyvalue(Sx.genset)) @ r(grass2)
    #             #coeff * r(base) @ (r(grass) * r(grass2))
    # for perm_pair in full_results:
    #     assert all(v > 0 for v in full_results[perm_pair].values()), f"Failure for perm_pair {perm_pair}, result {full_results[perm_pair]}"
    #     # assert results[perm_pair].change_basis(SchubertPolyBasis()).almosteq(Schub(perm_pair[0], num_parts) * Schub(perm_pair[1], num_parts)), f"Failure for perm_pair {perm_pair}, result {results[perm_pair]}"
    #     print(f"Success for perm_pair {perm_pair}, result {full_results[perm_pair]}")
    #     # key1 = KeyPoly.from_expr(sum([tt.factors[0].polyvalue(Sx.genset)*tt.factors[1].polyvalue(Sx.genset) for tt in t1.full_crystal]), length=num_parts)
    #     # key2 = KeyPoly.from_expr(sum([tt.factors[0].polyvalue(Sx.genset)*tt.factors[1].polyvalue(Sx.genset) for tt in t2.full_crystal]), length=num_parts)
    #     # assert all(v > 0 for v in key1.values()), f"Failure for hw {hw1}, tensor_hw {t1}, key {key1}"
    #     # assert all(v > 0 for v in key2.values()), f"Failure for hw {hw2}, tensor_hw {t2}, key {key2}"
    #     # key3 = key1 * key2
    #     # assert all(v > 0 for v in key3.values()), f"Failure for hw {hw1}, {hw2}, tensor_hw1 {t1}, tensor_hw2 {t2}, key1 {key1}, key2 {key2}, key3 {key3}"
    #     # hw_seen = {}
    #     # for base, factor_value in decomps.items():
    #     #     base_hw = base.to_highest_weight()[0]
    #     #     for other_rc in factor_value:
    #     #         tensor = CrystalGraphTensor(base, other_rc)
    #     #         hw = tensor.to_highest_weight()[0]
    #     #         if base_hw not in hw_seen:
    #     #             hw_seen[base_hw] = hw
    #     #         else:
    #     #             assert hw == hw_seen[base_hw], f"Failure for comp {comp}, base {base}, factor_value {factor_value}, decomps {decomps}"
            
    #         #assert all(v.almosteq(factor_value) for k, v in decomps.items() if k.perm == base.perm), f"Failure for comp {comp}, base {base}, factor_value {factor_value}, decomps {decomps}"
    #     print("Moraklistic success")

        