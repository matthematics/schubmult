from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *
from schubmult.utils.tuple_utils import pad_tuple
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
from sympy import pretty_print

def snap_key(rc):
    return next(iter([rcc for rcc in RCGraph.all_rc_graphs(rc.perm, len(rc), weight=rc.extremal_weight)]))

def forest_equal(rc1, rc2):
    return rc1.forest_weight == rc2.forest_weight and rc1.omega_invariant[1] == rc2.omega_invariant[1]


if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = ForestRCGraphRing()
    # Forest = PolynomialAlgebra(KeyPolyBasis(Sx.genset))
    # KeyDual = FreeAlgebra(KeyBasis)

    # for perm in perms:
    #     if perm.inv == 0:
    #         continue
    #     comp = pad_tuple(tuple(perm.trimcode), n - 1)
    #     forest_set = [rc for rc in RCGraph.all_rc_graphs(perm, n - 1) if rc.forest_weight == comp]

    #     for index in range(1, n - 1):
    #         fbranch = Forest(*comp).branch(index)
    #         tbranch = (Forest @ Forest).zero
    #         for rc in forest_set:
    #             rc1, rc2 = rc.vertical_cut(index)
    #             if rc1.forest_weight == rc1.length_vector and rc2.forest_weight == rc2.length_vector:
    #                 tbranch += Forest(*rc1.forest_weight) @ Forest(*rc2.forest_weight)
    #         assert fbranch.almosteq(tbranch), f"Failed for {perm} at index {index} with {len(forest_set)} RC graphs, got {tbranch} but expected {fbranch}"

    for perm1, perm2 in itertools.product(perms, repeat=2):
        comp1 = pad_tuple(tuple(perm1.trimcode), n - 1)
        comp2 = pad_tuple(tuple(perm2.trimcode), n - 1)

        forest_rc1 = r(RCGraph.principal_rc(perm1, n - 1))
        forest_rc2 = r(RCGraph.principal_rc(perm2, n - 1))
        # forest_rc1 = r.from_free_algebra_element(ForestDual(*comp1))
        # forest_rc2 = r.from_free_algebra_element(ForestDual(*comp2))
        # forest_rc1 = r(next(iter([rc for rc in RCGraph.all_rc_graphs(perm1, len(comp1)) if rc.forest_weight == comp1])))
        # forest_rc2 = r(next(iter([rc for rc in RCGraph.all_rc_graphs(perm2, len(comp2)) if rc.forest_weight == comp2])))

        rc_test_prod = ForestDual.zero
        #invariants_seen = {}
        rc_prod = (forest_rc1 * forest_rc2)
        # pretty_print(rc_prod)
        # print({k.forest_invariant: v for k, v in rc_prod.items()})
        for k, v in rc_prod.items():
            if k.forest_weight != pad_tuple(k.perm.trimcode, len(k)):
                raise ValueError(f"Failed for {comp1} * {comp2}, got non-forest RC graph {k} with weight {k.forest_weight} and perm {k.perm}, expected only forest RC graphs")
            rc_test_prod += v * ForestDual(*k.forest_weight)
        #rc_test_prd = sum([v * ForestDual(*k.forest_weight) for k, v in (r(forest_rc1) * r(forest_rc2)).items() if k.is_principal])
        #rc_test_prod = rc_prod.to_free_algebra_element(ForestBasis)
        expected_prd = ForestDual(*comp1) * ForestDual(*comp2)
        assert rc_test_prod.almosteq(expected_prd), f"Failed for {comp1} * {comp2}, got {rc_test_prod} but expected {expected_prd}"
        print(f"Passed for {comp1} * {comp2}")


        # key_set = [rc for rc in RCGraph.all_rc_graphs(perm, n - 1) if rc.extremal_weight == comp]

        # for index in range(1, n - 1):
        #     fbranch = Key(*comp).branch(index)
        #     tbranch = (Key @ Key).zero
        #     for rc in key_set:
        #         rc1, rc2 = rc.vertical_cut(index)
        #         if rc1.extremal_weight == rc1.length_vector and rc2.extremal_weight == rc2.length_vector:
        #             tbranch += Key(*rc1.extremal_weight) @ Key(*rc2.extremal_weight)
        #     assert fbranch.almosteq(tbranch), f"Failed for {perm} at index {index} with {len(key_set)} RC graphs, got {tbranch} but expected {fbranch}"
