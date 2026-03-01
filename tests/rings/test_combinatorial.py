def test_rc_graph_squash_product():
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

    ring = RCGraphRing()
    for max_d in range(n):
        build_product = {}
        for g_rc in grass_rcs[max_d]:
            for j in range(max_d + 1):
                for rc in rcs[j]:
                    rc2 = rc.resize(max_d)
                    the_rc = ring(rc2) % ring(g_rc)
                    for spank_rc, coeff in the_rc.items():
                        build_product[(g_rc.perm, rc2.perm)] = build_product.get((g_rc.perm, rc2.perm), ring.zero) + coeff * ring(spank_rc)

        for (g_perm, rc_perm), coeff in build_product.items():
            prod = Sx(g_perm) * Sx(rc_perm)
            assert all(prod.get(rc_result.perm, 0) == c for rc_result, c in coeff.items()), f"Error: Squash product mismatch for permutations {g_perm} and {rc_perm}, {coeff=} {prod=}"


def test_rc_graph_agrees_with_free_algebra():
    import itertools

    from schubmult import RCGraph, RCGraphRing, ASx, Permutation

    n = 4
    perms = Permutation.all_permutations(n)
    ring = RCGraphRing()

    for perm1, perm2 in itertools.product(perms, repeat=2):
        for len1 in range(len(perm1.trimcode), n):
            for len2 in range(len(perm2.trimcode), n):
                for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, len1), RCGraph.all_rc_graphs(perm2, len2)):
                    rc_elem = ring(rc1) * ring(rc2)
                    free_algebra_elem = ASx(rc1.perm, len(rc1)) * ASx(rc2.perm, len(rc2))
                    free_elem_test = rc_elem.to_free_algebra_element()
                    assert all(v == 0 for v in (free_algebra_elem - free_elem_test).values())


def test_rc_graph_coproduct():
    import itertools

    from schubmult import RCGraph, RCGraphRing, ASx, Permutation, SchubertBasis

    n = 4
    perms = Permutation.all_permutations(n)
    ring = RCGraphRing()

    for perm in perms:
        for length in range(len(perm.trimcode), n):
            for rc in RCGraph.all_rc_graphs(perm, length):
                rc_elem = ring(rc)
                free_algebra_elem = ASx(rc.perm, len(rc)).coproduct()
                rc_coprod = rc_elem.coproduct()
                assert all([v == 0 for v in (free_algebra_elem - sum([coeff * ASx(rc1.perm, len(rc1))@ASx(rc2.perm, len(rc2)) for (rc1, rc2), coeff in rc_coprod.items()])).values()])


def test_rc_bpd_ring_multiplication():
    import itertools

    from schubmult import RCGraph, BPD, Permutation, RCGraphRing, BPDRing
    from schubmult.symbolic import S

    n = 3
    perms = Permutation.all_permutations(n)

    ring1 = RCGraphRing()
    ring2 = BPDRing()
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for len1 in range(len(perm1.trimcode), n):
            for len2 in range(len(perm2.trimcode), n):
                for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, len1), RCGraph.all_rc_graphs(perm2, len2)):
                    rc_elem = ring1(rc1) * ring1(rc2)
                    bpd_elem = ring2(BPD.from_rc_graph(rc1)) * ring2(BPD.from_rc_graph(rc2))
                    assert all(v == S.Zero for v in (rc_elem - bpd_elem.to_rc_graph_ring_element()).values()), (
                        f"Error: RC graph ring element multiplication mismatch for permutations {perm1}, {perm2}:\n"
                        f"RC1:\n{rc1}\nRC2:\n{rc2}\nRC elem:\n{rc_elem}\nBPD elem:\n{bpd_elem.to_rc_graph_ring_element()}"
                    )


def test_bpd_agrees_with_free_algebra():
    import itertools

    from schubmult import BPD, BPDRing, ASx, Permutation

    n = 3
    perms = Permutation.all_permutations(n)
    ring = BPDRing()
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for len1 in range(len(perm1.trimcode), n):
            for len2 in range(len(perm2.trimcode), n):
                for rc1, rc2 in itertools.product(BPD.all_bpds(perm1, len1), BPD.all_bpds(perm2, len2)):
                    bpd_elem = ring(rc1) * ring(rc2)
                    free_algebra_elem = ASx(rc1.perm, len(rc1)) * ASx(rc2.perm, len(rc2))
                    free_elem_test = bpd_elem.to_free_algebra_element()
                    assert all(v == 0 for v in (free_algebra_elem - free_elem_test).values())


def test_eg_plactic_to_from_rc_graph():
    from schubmult import EGPlacticRing, RCGraph, RCGraphRing

    rc_ring = RCGraphRing()
    eg_ring = EGPlacticRing()
    rc = RCGraph([(4, 1), (3, 2), (), (4,)])
    rc_elem = rc_ring(rc)
    assert rc_elem.almosteq(eg_ring.from_rc_graph(rc).to_rc_graph_ring_element())


def test_eg_plactic_product_iso():
    from schubmult import EGPlacticRing, RCGraph, RCGraphRing

    rc_ring = RCGraphRing()
    eg_ring = EGPlacticRing()
    rc = RCGraph([(4, 1), (3, 2), (), (4,)])
    rc2 = RCGraph([(1,), (4, 3, 2), (4,), ()])
    rc_elem = rc_ring(rc) * rc_ring(rc2)
    assert rc_elem.almosteq((eg_ring.from_rc_graph(rc) * eg_ring.from_rc_graph(rc2)).to_rc_graph_ring_element())
