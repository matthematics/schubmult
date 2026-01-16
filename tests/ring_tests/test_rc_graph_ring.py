def test_squash_product():
    from schubmult import Permutation, RCGraph, RCGraphRing, Sx
    n = 4
    perms = Permutation.all_permutations(n)
    rcs = []
    grass_rcs = []
    for i in range(n):
        rcs.append(set())
        grass_rcs.append(set())

    for perm in perms:
        k = len(perm.trimcode)
        rc_set = RCGraph.all_rc_graphs(perm, k)
        rcs[k].update(rc_set)
        if len(perm.descents()) == 1:
            grass_rcs[k].update(rc_set)
    R = RCGraphRing()
    for max_d in range(n):
        build_product = {}
        for g_rc in grass_rcs[max_d]:
            
            for j in range(max_d + 1):
                for rc in rcs[j]:
                    rc2 = rc.resize(max_d)
                    the_rc = R(rc2) % R(g_rc)
                    for spank_rc, coeff in the_rc.items():
                        build_product[(g_rc.perm, rc2.perm)] = build_product.get((g_rc.perm, rc2.perm), R.zero) + coeff * R(spank_rc)
                
        for (g_perm, rc_perm), coeff in build_product.items():
            prod = Sx(g_perm) * Sx(rc_perm)
            assert all([prod.get(rc_result.perm, 0) == c for rc_result, c in coeff.items()]), f"Error: Squash product mismatch for permutations {g_perm} and {rc_perm}, {coeff=} {prod=}"
    
def test_agrees_with_free_algebra():
    from schubmult import RCGraph, RCGraphRing, ASx, Permutation
    import itertools
    n = 4
    
    perms = Permutation.all_permutations(n)

    for perm1, perm2 in itertools.product(perms, repeat=2):
        for len1 in range(len(perm1.trimcode), n):
            for len2 in range(len(perm2.trimcode), n):
                for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, len1), RCGraph.all_rc_graphs(perm2, len2)):
                    R = RCGraphRing()
                    rc_elem = R(rc1) * R(rc2)

                    free_algebra_elem = ASx(rc1.perm, len(rc1)) * ASx(rc2.perm, len(rc2))
                    free_elem_test = rc_elem.to_free_algebra_element()
                    assert all(v == 0 for v in (free_algebra_elem - free_elem_test).values())