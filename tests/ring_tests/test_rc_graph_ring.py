def test_squash_product():
    from schubmult import Permutation, RCGraph, RCGraphRing, Sx
    n = 5
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
                    for rc, coeff in the_rc.items():
                        build_product[(g_rc.perm, rc2.perm)] = build_product.get((g_rc.perm, rc2.perm), R.zero) + coeff * R(rc)
                
        for (g_perm, rc_perm), coeff in build_product.items():
            prod = Sx(g_perm) * Sx(rc_perm)
            assert all([prod.get(rc_result.perm, 0) == c for rc_result, c in coeff.items()])
    