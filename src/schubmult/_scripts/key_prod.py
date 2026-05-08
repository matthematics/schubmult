from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *
from schubmult.rings.combinatorial.hw_rc_ring import HWRCGraphRing
from sympy import pretty_print

# def snap_key(rc):
#     return next(iter([rcc for rcc in RCGraph.all_rc_graphs(rc.perm, len(rc), weight=rc.extremal_weight)]))

# def key_equal(rc1, rc2):
#     return rc1.extremal_weight == rc2.extremal_weight and rc1.omega_invariant[1] == rc2.omega_invariant[1]


if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    fr = HWRCGraphRing()
    KeyDual = FreeAlgebra(KeyBasis)


    for perm1, perm2 in itertools.product(perms, repeat=2):
        comp1 = perm1.trimcode
        comp2 = perm2.trimcode
        
        for length1, length2 in itertools.product(range(len(perm1.trimcode), n), range(len(perm2.trimcode), n)):
            
            comp1 = tuple(perm1.pad_code(length1))
            comp2 = tuple(perm2.pad_code(length2))
            
            rc_test_prod = 0
            key_rc1 =  sum([fr(potato) for potato in {bab.to_highest_weight()[0] for bab in RCGraph.all_key_rcs(comp1)}])
            key_rc2 = sum([fr(potato) for potato in {bab.to_highest_weight()[0] for bab in RCGraph.all_key_rcs(comp2)}])
            
            rc_prod = (key_rc1 * key_rc2)
            rc_test_prod2 = 0
            invariant_by_weight = {}
                    
            for k, v in rc_prod.items():
                rc_test_prod += v * ASx(k.perm,len(k)).change_basis(KeyBasis)#KeyDual(*k.extremal_weight)
                if k.extremal_weight not in invariant_by_weight:
                    invariant_by_weight[k.extremal_weight] = k.to_highest_weight()[0]
                elif invariant_by_weight[k.extremal_weight] != k.to_highest_weight()[0]:
                    continue
                rc_test_prod2 += v * KeyDual(*k.extremal_weight)
    
            expected_prd = KeyDual(*comp1) * KeyDual(*comp2)
            assert rc_test_prod.almosteq(expected_prd), f"Failed for {comp1} * {comp2}, got {rc_test_prod} but expected {expected_prd}\n{rc_test_prod-expected_prd=}"
            assert rc_test_prod2.almosteq(expected_prd), f"Failed for {comp1} * {comp2}, got {rc_test_prod2} but expected {expected_prd}\n{rc_test_prod2-expected_prd=}"
            print(f"Passed for {comp1} * {comp2}")

    