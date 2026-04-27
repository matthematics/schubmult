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
    #     forest_set = [rc for rc in RCGraph.all_rc_graphs(perm, n - 1) if rc.forest_weight == comp]

    #     for index in range(1, n - 1):
    #         fbranch = Forest(*comp).branch(index)
    #         tbranch = (Forest @ Forest).zero
    #         for rc in forest_set:
    #             rc1, rc2 = rc.vertical_cut(index)
    #             if rc1.forest_weight == rc1.length_vector and rc2.forest_weight == rc2.length_vector:
    #                 tbranch += Forest(*rc1.forest_weight) @ Forest(*rc2.forest_weight)
    #         assert fbranch.almosteq(tbranch), f"Failed for {perm} at index {index} with {len(forest_set)} RC graphs, got {tbranch} but expected {fbranch}"
    #         print(f"Success {comp}")

    for perm1, perm2 in itertools.product(perms, repeat=2):
        for length1, length2 in itertools.product(range(len(perm1.trimcode), n), range(len(perm2.trimcode), n)):
            comp1 = perm1.pad_code(length1)
            comp2 = perm2.pad_code(length2)
            for forest_rc1_base, forest_rc2_base in itertools.product(RCGraph.all_forest_rcs(comp1, weight=comp1), RCGraph.all_forest_rcs(comp2, weight=comp2)):
                forest_rc1 = r(forest_rc1_base)
                forest_rc2 = r(forest_rc2_base)

                rc_test_prod = ForestDual.zero

                rc_prod = (forest_rc1 * forest_rc2)
                seen_by_perm = {}
                for k, v in rc_prod.items():
                    if k.forest_weight not in seen_by_perm:
                        seen_by_perm[k.forest_weight] = k.forest_invariant
                    elif seen_by_perm[k.forest_weight] != k.forest_invariant:
                        continue
                    # if k.forest_weight != k.perm.pad_code(len(k)):
                    #     continue
                    #     raise ValueError(f"Failed for {comp1} * {comp2}, got non-forest RC graph {k} with weight {k.forest_weight} and perm {k.perm}, expected only forest RC graphs")
                    rc_test_prod += v * ForestDual(*k.forest_weight)
                expected_prd = ForestDual(*comp1) * ForestDual(*comp2)
                assert rc_test_prod.almosteq(expected_prd), f"Failed for {comp1} * {comp2}, got {rc_test_prod} but expected {expected_prd}\n{rc_prod=}"
                print(f"Passed for {comp1} * {comp2}")
