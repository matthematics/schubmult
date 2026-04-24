from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import KeyBasis, FreeAlgebra
from schubmult.utils.tuple_utils import pad_tuple
from schubmult.rings.combinatorial.key_rc_ring import KeyRCGraphRing
from sympy import pretty_print

#KeyDual = FreeAlgebra(KeyBasis)
def snap_key(rc):
    return next(iter([rcc for rcc in RCGraph.all_rc_graphs(rc.perm, len(rc), weight=rc.extremal_weight)]))

# def key_equal(rc1, rc2):
#     return rc1.extremal_weight == rc2.extremal_weight and rc1.omega_invariant[1] == rc2.omega_invariant[1]


if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = KeyRCGraphRing()
    # Key = PolynomialAlgebra(KeyPolyBasis(Sx.genset))
    KeyDual = FreeAlgebra(KeyBasis)

    # for perm in perms:
    #     if perm.inv == 0:
    #         continue
    #     comp = pad_tuple(tuple(perm.trimcode), n - 1)
    #     key_set = [rc for rc in RCGraph.all_rc_graphs(perm, n - 1) if rc.extremal_weight == comp]

    #     for index in range(1, n - 1):
    #         fbranch = Key(*comp).branch(index)
    #         tbranch = (Key @ Key).zero
    #         for rc in key_set:
    #             rc1, rc2 = rc.vertical_cut(index)
    #             if rc1.extremal_weight == rc1.length_vector and rc2.extremal_weight == rc2.length_vector:
    #                 tbranch += Key(*rc1.extremal_weight) @ Key(*rc2.extremal_weight)
    #         assert fbranch.almosteq(tbranch), f"Failed for {perm} at index {index} with {len(key_set)} RC graphs, got {tbranch} but expected {fbranch}"

    for perm1, perm2 in itertools.product(perms, repeat=2):
        comp1 = pad_tuple(tuple(perm1.trimcode), n - 1)
        comp2 = pad_tuple(tuple(perm2.trimcode), n - 1)

        key_rc1 = r(RCGraph.principal_rc(perm1, n - 1))
        key_rc2 = r(RCGraph.principal_rc(perm2, n - 1))
        # key_rc1 = r.from_free_algebra_element(KeyDual(*comp1))
        # key_rc2 = r.from_free_algebra_element(KeyDual(*comp2))
        # key_rc1 = r(next(iter([rc for rc in RCGraph.all_rc_graphs(perm1, len(comp1)) if rc.extremal_weight == comp1])))
        # key_rc2 = r(next(iter([rc for rc in RCGraph.all_rc_graphs(perm2, len(comp2)) if rc.extremal_weight == comp2])))

        rc_test_prod = KeyDual.zero
        #invariants_seen = {}
        rc_prod = (key_rc1 * key_rc2)
        # pretty_print(rc_prod)
        # print({k.key_invariant: v for k, v in rc_prod.items()})
        for k, v in rc_prod.items():
            # if k.extremal_weight != pad_tuple(k.perm.trimcode, len(k)):
            #     raise ValueError(f"Failed for {comp1} * {comp2}, got non-forest RC graph {k} with weight {k.extremal_weight} and perm {k.perm}, expected only forest RC graphs")
            rc_test_prod += v * KeyDual(*k.extremal_weight)
        #rc_test_prd = sum([v * KeyDual(*k.extremal_weight) for k, v in (r(key_rc1) * r(key_rc2)).items() if k.is_principal])
        #rc_test_prod = rc_prod.to_free_algebra_element(KeyBasis)
        expected_prd = KeyDual(*comp1) * KeyDual(*comp2)
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
