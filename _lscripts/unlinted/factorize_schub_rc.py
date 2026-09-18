from schubmult import *
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
    grass_hw_seen = {}
    for comp in comps:
        print(f"Testing comp {comp}...")
        perm = uncode(comp)
        if len(perm) < num_parts + 1:
            continue
        #decomps = {}
        
        for rc in RCGraph.all_rc_graphs(perm, num_parts):
            base, grass =  rc.squash_decomp()
            tensor = CrystalGraphTensor(base, grass)
            tensor_hw = tensor.to_highest_weight()[0]
            hw = rc.to_highest_weight()[0]
            if hw not in hw_seen:
                hw_seen[hw] = tensor_hw
            else:
                assert tensor_hw == hw_seen[hw], f"Failure for comp {comp}, base {base}, grass {grass}, hw {hw}, tensor_hw {tensor_hw}"
            if perm.inv == 0 or perm.descents() == {num_parts}:
                if hw not in grass_hw_seen:
                    grass_hw_seen[hw] = tensor_hw
                
            #decomps[base] = decomps.get(base, r.zero) + r(grass)
    results = {}
    for hw1, hw2 in itertools.product(hw_seen.keys(), grass_hw_seen.keys()):
        t1 = hw_seen[hw1]
        t2 = hw_seen[hw2]
        results[(hw1, hw2)] = results.get((hw1, hw2), 0)
        for tt1 in t1.full_crystal:
            for tt2 in t2.full_crystal:
                # we have tt1.factors[0] @ tt1.factors[1] * tt2.factors[0] @ tt2.factors[1]
                middle = tt1.factors[1].left_squash(tt2.factors[0])
                base, grass2 = middle.squash_decomp()
                right = grass2.squash_product(tt2.factors[1])
                # #results[(hw1.perm, hw2.perm)] = results.get((hw1.perm, hw2.perm), 0) + KeyPoly.from_expr(tt1.factors[0].polyvalue(Sx.genset) * tt2.factors[0].polyvalue(Sx.genset), length=num_parts) @ (r(tt1.factors[1])*r(tt2.factors[1]))
                
                # #results[(hw1, hw2)] += r(tt1.factors[0]) @ r(base) @ right
                
                results[(hw1, hw2)] += r(tt1.factors[0]) @ (r(base) * r(right))
                #Schub.from_expr(tt1.factors[0].polyvalue(Sx.genset) * tt2.factors[0].polyvalue(Sx.genset), length=num_parts) @ (r(tt1.factors[1])*r(tt2.factors[1]))
                #+ KeyPoly.from_expr(tt1.factors[0].polyvalue(Sx.genset) * tt1.factors[1].polyvalue(Sx.genset), length=num_parts)@ KeyPoly.from_expr((r(tt2.factors[0])*r(tt2.factors[1])).polyvalue(Sx.genset), length=num_parts)
    int_results = {}
    for (hw1, hw2), result in results.items():
        for (key1, key2, key3), coeff in result.items():
            
            #int_results[(hw1.perm, hw2.perm)] = int_results.get((hw1.perm, hw2.perm), 0) + Schub.from_expr(key1.polyvalue(Sx.genset) * key2.polyvalue(Sx.genset), length=num_parts) @ r(key3)
            #int_results[(hw1.perm, hw2.perm)] = int_results.get((hw1.perm, hw2.perm), 0) + coeff * KeyPoly.from_expr(key1.polyvalue(Sx.genset)key3.polyvalue(Sx.genset), length=num_parts)
            int_results[(hw1.perm, hw2.perm)] = int_results.get((hw1.perm, hw2.perm), 0) + coeff * Schub.from_expr(key1.polyvalue(Sx.genset)*key2.polyvalue(Sx.genset), length=num_parts) @ r(key3)
            #, length=num_parts) @ KeyPoly.from_expr(key3.polyvalue(Sx.genset), length=num_parts)
    #full_results = {
    full_results = {}
    for (perm1, perm2), result in int_results.items():
        for ((perm, length), grass_rc), coeff in result.items():
            for new_rc in RCGraph.all_rc_graphs(perm, num_parts):
                base, grass2_rc = new_rc.squash_decomp()
                #fabrika = r(base) * r(grass2_rc) * r(grass_rc)
                full_results[(perm1, perm2)] = full_results.get((perm1, perm2), 0) + coeff * KeyPoly.from_expr(base.polyvalue(Sx.genset)) @ (r(grass2_rc) * r(grass_rc))
                    #base) @ (r(grass2_rc) * r(grass_rc))
        # for (base1, base2, grass2), coeff in result.items():
        #     # for new_rc in RCGraph.all_rc_graphs(perm, num_parts):
        #     #     base, grass = new_rc.squash_decomp()
        #     full_results[(hw1, hw2)] = full_results.get((hw1, hw2), 0) + coeff * Schub.from_expr(base1.polyvalue(Sx.genset) * base2.polyvalue(Sx.genset)) @ r(grass2)
                #coeff * r(base) @ (r(grass) * r(grass2))
    for perm_pair in full_results:
        assert all(v > 0 for v in full_results[perm_pair].values()), f"Failure for perm_pair {perm_pair}, result {full_results[perm_pair]}"
        # assert results[perm_pair].change_basis(SchubertPolyBasis()).almosteq(Schub(perm_pair[0], num_parts) * Schub(perm_pair[1], num_parts)), f"Failure for perm_pair {perm_pair}, result {results[perm_pair]}"
        print(f"Success for perm_pair {perm_pair}, result {full_results[perm_pair]}")
        # key1 = KeyPoly.from_expr(sum([tt.factors[0].polyvalue(Sx.genset)*tt.factors[1].polyvalue(Sx.genset) for tt in t1.full_crystal]), length=num_parts)
        # key2 = KeyPoly.from_expr(sum([tt.factors[0].polyvalue(Sx.genset)*tt.factors[1].polyvalue(Sx.genset) for tt in t2.full_crystal]), length=num_parts)
        # assert all(v > 0 for v in key1.values()), f"Failure for hw {hw1}, tensor_hw {t1}, key {key1}"
        # assert all(v > 0 for v in key2.values()), f"Failure for hw {hw2}, tensor_hw {t2}, key {key2}"
        # key3 = key1 * key2
        # assert all(v > 0 for v in key3.values()), f"Failure for hw {hw1}, {hw2}, tensor_hw1 {t1}, tensor_hw2 {t2}, key1 {key1}, key2 {key2}, key3 {key3}"
        # hw_seen = {}
        # for base, factor_value in decomps.items():
        #     base_hw = base.to_highest_weight()[0]
        #     for other_rc in factor_value:
        #         tensor = CrystalGraphTensor(base, other_rc)
        #         hw = tensor.to_highest_weight()[0]
        #         if base_hw not in hw_seen:
        #             hw_seen[base_hw] = hw
        #         else:
        #             assert hw == hw_seen[base_hw], f"Failure for comp {comp}, base {base}, factor_value {factor_value}, decomps {decomps}"
            
            #assert all(v.almosteq(factor_value) for k, v in decomps.items() if k.perm == base.perm), f"Failure for comp {comp}, base {base}, factor_value {factor_value}, decomps {decomps}"
        print("Moraklistic success")

        