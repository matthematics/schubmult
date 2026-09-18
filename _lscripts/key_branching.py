from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import KeyBasis, FreeAlgebra
from schubmult.utils.tuple_utils import pad_tuple
from schubmult.rings.combinatorial.hw_rc_ring import HWRCGraphRing
from sympy import pretty_print

#KeyDual = FreeAlgebra(KeyBasis)
# def snap_key(rc):
#     return next(iter([rcc for rcc in RCGraph.all_rc_graphs(rc.perm, len(rc), weight=rc.extremal_weight)]))

# def key_equal(rc1, rc2):
#     return rc1.extremal_weight == rc2.extremal_weight and rc1.omega_invariant[1] == rc2.omega_invariant[1]


if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = HWRCGraphRing()
    # r = RCGraphRing()
    Key = PolynomialAlgebra(KeyPolyBasis(Sx.genset))
    KeyDual = FreeAlgebra(KeyBasis)

    # for perm in perms:
    #     if perm.inv == 0:
    #         continue
    #     comp = pad_tuple(tuple(perm.trimcode), n - 1)
    #     # key_set = [rc for rc in RCGraph.all_rc_graphs(perm, n - 1) if rc.extremal_weight == comp]
    #     key_set = RCGraph.principal_rc(perm, n-1).full_crystal
    #     for index in range(1, n - 1):
    #         fbranch = Key(*comp).branch(index)
    #         tbranch = (Key @ Key).zero
    #         for rc in key_set:
    #             rc1, rc2 = rc.vertical_cut(index)
    #             if rc1.extremal_weight == rc1.length_vector and rc2.extremal_weight == rc2.length_vector:
    #                 tbranch += Key(*rc1.extremal_weight) @ Key(*rc2.extremal_weight)
    #         assert fbranch.almosteq(tbranch), f"Failed for {perm} at index {index} with {len(key_set)} RC graphs, got {tbranch} but expected {fbranch}"
    #     print(f"Success {comp}")

    for perm1, perm2 in itertools.product(perms, repeat=2):
        for length1, length2 in itertools.product(range(max(1,len(perm1.trimcode)), n), range(max(1,len(perm2.trimcode)), n)):
            comp1 = pad_tuple(tuple(perm1.trimcode), length1)
            comp2 = pad_tuple(tuple(perm2.trimcode), length2)
            for key_rc1_base, key_rc2_base in itertools.product(RCGraph.all_hw_rcs(uncode(comp1), length1), RCGraph.all_hw_rcs(uncode(comp2), length2)):
            #if True:
                # key_rc1_base = RCGraph.principal_rc(perm1, length1)    
                # key_rc2_base = RCGraph.principal_rc(perm2, length2)    
                if key_rc1_base.extremal_weight != comp1 or key_rc2_base.extremal_weight != comp2:
                    continue
                key_rc1 = r(key_rc1_base)
                key_rc2 = r(key_rc2_base)

                rc_test_prod = KeyDual.zero

                rc_prod = (key_rc1 * key_rc2)
                extw_seen = {}
                for k, v in rc_prod.items():
                    # if k.extremal_weight != k.perm.pad_code(len(k)):
                    #     continue
                    #     raise ValueError(f"Failed for {comp1} * {comp2}, got non-forest RC graph {k} with weight {k.extremal_weight} and perm {k.perm}, expected only forest RC graphs")
                    # if k.extremal_weight == k.length_vector:
                    # if k.extremal_weight not in extw_seen:
                    #     extw_seen[k.extremal_weight] = k.to_highest_weight()[0]
                    # elif extw_seen[k.extremal_weight] != k.to_highest_weight()[0]:
                    #     continue
                    rc_test_prod += v * KeyDual(*k.extremal_weight)
                    #    hw_seen.add(k.to_highest_weight()[0])
                expected_prd = KeyDual(*comp1) * KeyDual(*comp2)
                assert rc_test_prod.almosteq(expected_prd), f"Failed for {comp1} * {comp2}, got {rc_test_prod} but expected {expected_prd}\n{rc_prod=}"
                print(f"Passed for {comp1} * {comp2}")
            

    # #     # key_set = [rc for rc in RCGraph.all_rc_graphs(perm, n - 1) if rc.extremal_weight == comp]

    # #     # for index in range(1, n - 1):
    # #     #     fbranch = Key(*comp).branch(index)
    # #     #     tbranch = (Key @ Key).zero
    # #     #     for rc in key_set:
    # #     #         rc1, rc2 = rc.vertical_cut(index)
    # #     #         if rc1.extremal_weight == rc1.length_vector and rc2.extremal_weight == rc2.length_vector:
    # #     #             tbranch += Key(*rc1.extremal_weight) @ Key(*rc2.extremal_weight)
    # #     #     assert fbranch.almosteq(tbranch), f"Failed for {perm} at index {index} with {len(key_set)} RC graphs, got {tbranch} but expected {fbranch}"
