from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
from sympy import pretty_print

# def snap_key(rc):
#     return next(iter([rcc for rcc in RCGraph.all_rc_graphs(rc.perm, len(rc), weight=rc.extremal_weight)]))

# def forest_equal(rc1, rc2):
#     return rc1.forest_weight == rc2.forest_weight and rc1.omega_invariant[1] == rc2.omega_invariant[1]


if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    # r = ForestRCGraphRing()
    r = RCGraphRing()
    Forest = PolynomialAlgebra(ForestPolyBasis(Sx.genset))
    ForestDual = FreeAlgebra(ForestBasis)



    # for perm in perms:
    #     if perm.inv == 0:
    #         continue
    #     comp = perm.pad_code(n-1)
    #     forest_set = [rc for rc in RCGraph.all_rc_graphs(perm, n - 1)]
        
    #     for index in range(1, n - 1):
    #         fbranch = Forest(*comp).branch(index)
    #         #tbranch = (Forest @ Forest).zero
    #         by_invariant = {}
    #         for rc in forest_set:
    #             rc1, rc2 = rc.vertical_cut(index)
    #             if rc1.forest_weight == rc1.length_vector and rc2.forest_weight == rc2.length_vector:
    #                 by_invariant[rc.forest_invariant] = by_invariant.get(rc.forest_invariant, 0) + Forest(*rc1.forest_weight) @ Forest(*rc2.forest_weight)
    #         for invar, val in by_invariant.items():
    #             assert fbranch.almosteq(val), f"Failed for {perm} at index {index} with {len(forest_set)} RC graphs, got {val} but expected {fbranch}"
    #         print(f"Success {comp}")

    for perm1, perm2 in itertools.product(perms, repeat=2):
        if perm1.inv == 0 and perm2.inv == 0:
            continue
        comp1 = perm1.trimcode
        comp2 = perm2.trimcode
        # forest_rc1_base = r.from_dict({k: 1 for k, v in r.from_free_algebra_element(ForestDual(*comp1)).items() if k.forest_weight == tuple(comp1) and (k.perm, len(k)) in ForestDual(*comp1).change_basis(SchubertBasis) and v == 1})
        
        # forest_rc2_base = r.from_dict({k: 1 for k, v in r.from_free_algebra_element(ForestDual(*comp2)).items() if k.forest_weight == tuple(comp2) and (k.perm, len(k)) in ForestDual(*comp2).change_basis(SchubertBasis) and v == 1})
        # print(len(forest_rc1_base))
        # print(len(forest_rc2_base))
        for length1, length2 in itertools.product(range(max(1,len(perm1.trimcode)), n), range(max(1,len(perm2.trimcode)), n)):
            
            comp1 = perm1.pad_code(length1)
            comp2 = perm2.pad_code(length2)
            
            
            
            for rc11 in RCGraph.all_forest_rcs(comp1):
                forest_rc1 = r(rc11)
                for rc22 in RCGraph.all_forest_rcs(comp2):
                    rc_test_prod = ForestDual.zero    
                    # forest_rc1 = forest_rc1_base.resize(length1)
                    # forest_rc2 = forest_rc2_base.resize(length2)
                    
                    
                    forest_rc2 = r(rc22)
                
                # if rc11.forest_weight not in seen_by_weight:
                #     seen_by_weight[rc11.forest_weight] = rc11.forest_invariant
                # elif seen_by_weight[rc11.forest_weight] != rc11.forest_invariant:
                #     continue
                # if rc22.forest_weight not in seen_by_weight:
                #     seen_by_weight[rc22.forest_weight] = rc22.forest_invariant
                # elif seen_by_weight[rc22.forest_weight] != rc22.forest_invariant:
                #     continue
                    #seen_by_weight = {}
                    seen = set()
                    # forest1 = ForestDual(*comp1).change_basis(SchubertBasis)
                    # forest2 = ForestDual(*comp2).change_basis(SchubertBasis)
                    # #for forest_rc1_base, forest_rc2_base in itertools.product(RCGraph.all_forest_rcs(comp1), RCGraph.all_forest_rcs(comp2)):
                    # forest_rc1 = sum([r(RCGraph.principal_rc(perm, length)) for (perm, length), v in forest1.items() if v != 0])
                    # forest_rc2 = sum([r(RCGraph.principal_rc(perm, length)) for (perm, length), v in forest2.items() if v != 0])
                    
                    # assert (forest_rc1_base.perm, len(forest_rc1_base)) in ForestDual(*comp1).change_basis(SchubertBasis), f"Failed for {comp1}, got non-forest RC graph {forest_rc1} with weight {forest_rc1.forest_weight} and perm {forest_rc1.perm}, expected only forest RC graphs"
                    # assert (forest_rc2_base.perm, len(forest_rc2_base)) in ForestDual(*comp2).change_basis(SchubertBasis), f"Failed for {comp2}, got non-forest RC graph {forest_rc2} with weight {forest_rc2.forest_weight} and perm {forest_rc2.perm}, expected only forest RC graphs"
                    

                    rc_prod = (forest_rc1 * forest_rc2)
                    
                    # for rc1 in forest_rc1:
                    #     seen_by_weight[rc1.forest_invariant] = rc1
                    # for rc2 in forest_rc2:
                    #     #if rc2.forest_invariant in seen_by_weight and 
                    #     seen_by_weight[rc2.forest_invariant] = rc2
                    for k, v in rc_prod.items():
                        # if k.forest_weight not in seen_by_weight:
                        #     seen_by_weight[k.forest_weight] = k.forest_invariant
                        # elif seen_by_weight[k.forest_weight] != k.forest_invariant:
                        #     continue
                        # if k.forest_weight != k.perm.pad_code(len(k)):
                        #     continue
                        #     raise ValueError(f"Failed for {comp1} * {comp2}, got non-forest RC graph {k} with weight {k.forest_weight} and perm {k.perm}, expected only forest RC graphs")
                        #if (k.perm, len(k)) in ForestDual(*k.forest_weight).change_basis(SchubertBasis):
                        # if k.forest_weight not in seen_by_weight:
                        #     # if k.forest_weight == comp1 and k != next(iter(forest_rc1_base)):
                        #     #     continue
                        #     # if k.forest_weight == comp2 and k != next(iter(forest_rc2_base)):
                        #     #     continue
                        #     seen_by_weight[k.forest_weight] = k.forest_invariant
                        # elif seen_by_weight[k.forest_weight] != k.forest_invariant:
                        #     continue
                        if k.forest_weight in seen:#invariant in seen:
                            continue
                        seen.add(k.forest_weight)#invariant)
                        #rc_test_prod += v * ASx(k.perm,len(k)).change_basis(ForestBasis)#ForestDual(*k.forest_weight)
                        rc_test_prod += v * ForestDual(*k.forest_weight)
                        
                expected_prd = ForestDual(*comp1) * ForestDual(*comp2)
                assert rc_test_prod.almosteq(expected_prd), f"Failed for {comp1} * {comp2}, got {rc_test_prod} but expected {expected_prd}\n{rc_test_prod-expected_prd=}"
            print(f"Passed for {comp1} * {comp2}")
