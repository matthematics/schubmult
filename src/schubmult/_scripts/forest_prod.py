from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing, _canonical_rc
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
    fr = ForestRCGraphRing()
    
    ForestDual = FreeAlgebra(ForestBasis)

    for perm1, perm2 in itertools.product(perms, repeat=2):
        comp1 = perm1.trimcode
        comp2 = perm2.trimcode
        
        for length1, length2 in itertools.product(range(len(perm1.trimcode), n), range(len(perm2.trimcode), n)):
            
            comp1 = tuple(perm1.pad_code(length1))
            comp2 = tuple(perm2.pad_code(length2))
            
            rc_test_prod = 0
            forest_rc1 =  sum([fr(potato) for potato in {_canonical_rc(bab) for bab in RCGraph.all_forest_rcs(comp1)}])
            forest_rc2 = sum([fr(potato) for potato in {_canonical_rc(bab) for bab in RCGraph.all_forest_rcs(comp2)}])
            
            rc_prod = (forest_rc1 * forest_rc2)
            rc_test_prod2 = 0
            invariant_by_weight = {}
                    
            for k, v in rc_prod.items():
                rc_test_prod += v * ASx(k.perm,len(k)).change_basis(ForestBasis)#ForestDual(*k.forest_weight)
                if k.forest_weight not in invariant_by_weight:
                    invariant_by_weight[k.forest_weight] = k.forest_invariant
                elif invariant_by_weight[k.forest_weight] != k.forest_invariant:
                    continue
                rc_test_prod2 += v * ForestDual(*k.forest_weight)
    
            expected_prd = ForestDual(*comp1) * ForestDual(*comp2)
            assert rc_test_prod.almosteq(expected_prd), f"Failed for {comp1} * {comp2}, got {rc_test_prod} but expected {expected_prd}\n{rc_test_prod-expected_prd=}"
            assert rc_test_prod2.almosteq(expected_prd), f"Failed for {comp1} * {comp2}, got {rc_test_prod2} but expected {expected_prd}\n{rc_test_prod2-expected_prd=}"
            print(f"Passed for {comp1} * {comp2}")
