import pytest
def test_rc_graph_squash_product():
    from schubmult import Permutation, RCGraph, DualRCGraphRing, Sx

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

    ring = DualRCGraphRing()
    
    for max_d in range(n):
        build_product = {}
        bad_perm_pairs = set()
        for g_rc in grass_rcs[max_d]:
            for j in range(max_d + 1):
                for rc in rcs[j]:
                    rc2 = rc.resize(max_d)
                    try:
                        the_rc = ring(rc2.squash_product(g_rc))
                    except NotImplementedError:
                        bad_perm_pairs.add((g_rc.perm, rc2.perm))
                        continue
                    for spank_rc, coeff in the_rc.items():
                        build_product[(g_rc.perm, rc2.perm)] = build_product.get((g_rc.perm, rc2.perm), ring.zero) + coeff * ring(spank_rc)

        for (perm1, perm2) in bad_perm_pairs:
            if (perm1, perm2) in build_product:
                del build_product[(perm1, perm2)]

        for (g_perm, rc_perm), coeff in build_product.items():
            # if (g_perm, rc_perm) in bad_perm_pairs or (rc_perm, g_perm) in bad_perm_pairs:
            #     continue
            prod = Sx(g_perm) * Sx(rc_perm)
            assert all(prod.get(rc_result.perm, 0) == c for rc_result, c in coeff.items()), f"Error: Squash product mismatch for permutations {g_perm} and {rc_perm}, {coeff=} {prod=}"

def test_rc_graph_left_squash_product():
    from schubmult import Permutation, RCGraph, DualRCGraphRing, Sx

    ring = DualRCGraphRing()

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

    
    for max_d in range(n):
        build_product = {}
        bad_perm_pairs = set()
        for g_rc in grass_rcs[max_d]:
            for j in range(max_d + 1):
                for rc in rcs[j]:
                    rc2 = rc.resize(max_d)
                    try:
                        the_rc = rc2.squash_product(g_rc)
                    except NotImplementedError:
                        bad_perm_pairs.add((g_rc.perm, rc2.perm))
                        continue
                    build_product[(g_rc.perm, rc2.perm)] = build_product.get((g_rc.perm, rc2.perm), ring.zero) + ring(the_rc)

        for (perm1, perm2) in bad_perm_pairs:
            if (perm1, perm2) in build_product:
                del build_product[(perm1, perm2)]

        for (g_perm, rc_perm), coeff in build_product.items():
            # if (g_perm, rc_perm) in bad_perm_pairs or (rc_perm, g_perm) in bad_perm_pairs:
            #     continue
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


def test_rc_graph_empty_list_constructor_has_zero_rows():
    from schubmult import RCGraph

    assert len(RCGraph([])) == 0


def test_rc_graph_ring_resize_resizes_empty_graph_too():
    from schubmult import RCGraph
    from schubmult.rings.combinatorial.rc_graph_ring import RCGraphRing

    ring = RCGraphRing()
    resized = ring(RCGraph([])).resize(3)

    assert len(resized) == 1
    only_key = next(iter(resized.keys()))
    assert len(only_key) == 3
    assert resized[only_key] == 1


def test_rc_graph_ring_empty_basis_not_scalar_one_term():
    from schubmult import RCGraph
    from schubmult.rings.combinatorial.rc_graph_ring import RCGraphRing
    from schubmult.symbolic import S

    ring = RCGraphRing()
    term = ring(RCGraph([])).as_ordered_terms()[0]

    assert term != S.One


def test_forest_rc_graph_ring_empty_basis_not_scalar_one_term():
    from schubmult import RCGraph
    from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
    from schubmult.symbolic import S

    ring = ForestRCGraphRing()
    term = ring(RCGraph([])).as_ordered_terms()[0]

    assert term != S.One


def test_tensor_ring_empty_basis_not_scalar_one_term():
    from schubmult import RCGraph
    from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
    from schubmult.symbolic import S

    ring = ForestRCGraphRing() @ ForestRCGraphRing()
    term = ring((RCGraph([]), RCGraph([]))).as_ordered_terms()[0]

    assert term != S.One


def test_wc_graph_ring_product_agrees_with_free_algebra():
    import itertools

    from schubmult import AGx, Permutation
    from schubmult.combinatorics.wc_graph import WCGraph
    from schubmult.rings.combinatorial.wc_graph_ring import WCGraphRing

    def to_agx(elem):
        res = AGx.zero
        for wc, coeff in elem.items():
            res += coeff * AGx(wc.perm, len(wc))
        return res

    n = 3
    perms = Permutation.all_permutations(n)
    ring = WCGraphRing()

    for perm1, perm2 in itertools.product(perms, repeat=2):
        for len1 in range(len(perm1.trimcode), n):
            for len2 in range(len(perm2.trimcode), n):
                for wc1, wc2 in itertools.product(WCGraph.all_wc_graphs(perm1, len1), WCGraph.all_wc_graphs(perm2, len2)):
                    wc_elem = ring(wc1) * ring(wc2)
                    free_algebra_elem = AGx(wc1.perm, len(wc1)) * AGx(wc2.perm, len(wc2))
                    free_elem_test = to_agx(wc_elem)
                    assert all(v == 0 for v in (free_algebra_elem - free_elem_test).values()), (
                        f"WCGraphRing product disagrees with AGx for {perm1}, {perm2}: diff={free_algebra_elem - free_elem_test} {wc_elem=} {free_algebra_elem=} {free_elem_test=}"
                    )


def test_wc_graph_ring_groth_polyvalue_is_grothendieck_poly():
    import itertools

    from schubmult import Permutation
    from schubmult.abc import x
    from schubmult.rings.combinatorial.wc_graph_ring import WCGraphRing
    from schubmult.rings.schubert.grothendieck_ring import Gx
    from schubmult.symbolic import S, expand
    from schubmult.symbolic.common_polys import grothendieck_poly
    from schubmult.symbolic.poly.variables import ZeroGeneratingSet

    ring = WCGraphRing()
    beta = Gx._beta
    zz = ZeroGeneratingSet()

    for pl in itertools.permutations(range(1, 4)):
        w = Permutation(list(pl))
        elem = ring.groth(w)
        polyval = S.Zero
        for wc, coeff in elem.items():
            polyval += coeff * wc.polyvalue(x, beta=beta, prop_beta=True)
        direct = grothendieck_poly(w, x, zz, beta)
        assert expand(polyval - direct, deep=True) == S.Zero, f"Groth polyvalue mismatch for {w}"


def test_wc_graph_ring_one_is_empty_graph():
    from schubmult.combinatorics.wc_graph import WCGraph
    from schubmult.rings.combinatorial.wc_graph_ring import WCGraphRing

    ring = WCGraphRing()
    keys = list(ring.one.keys())

    assert keys == [WCGraph([])]
    assert ring.one[WCGraph([])] == 1


def test_wc_graph_ring_empty_basis_not_scalar_one_term():
    from schubmult.combinatorics.wc_graph import WCGraph
    from schubmult.rings.combinatorial.wc_graph_ring import WCGraphRing
    from schubmult.symbolic import S

    ring = WCGraphRing()
    term = ring(WCGraph([])).as_ordered_terms()[0]

    assert term != S.One


def test_wc_graph_ring_resize_resizes_empty_graph_too():
    from schubmult.combinatorics.wc_graph import WCGraph
    from schubmult.rings.combinatorial.wc_graph_ring import WCGraphRing

    ring = WCGraphRing()
    resized = ring(WCGraph([])).resize(3)

    assert len(resized) == 1
    only_key = next(iter(resized.keys()))
    assert len(only_key) == 3
    assert resized[only_key] == 1

def test_bounded_rc_factor_algebra_multiplies_schuberts():
    from schubmult import Permutation, BoundedRCFactorAlgebra, Sx

    n = 3
    perms = Permutation.all_permutations(n)
    ring = BoundedRCFactorAlgebra()

    for perm1 in perms:
        for perm2 in perms:
            elem1 = ring.full_schub_elem(perm1, n + 2)
            elem2 = ring.full_schub_elem(perm2, n + 2)
            prod = (elem1 * elem2).to_rc_graph_ring_element()
            expected_prod = Sx(perm1) * Sx(perm2)
            assert all(expected_prod.get(rc_result.perm, 0) == c for rc_result, c in prod.items()), f"Error: Bounded RC factor algebra multiplication mismatch for permutations {perm1} and {perm2}, {prod=} {expected_prod=}"

def test_bounded_wc_factor_algebra_multiplies_grothendiecks():
    from schubmult import Permutation, BoundedWCFactorAlgebra, Gx

    n = 3
    perms = Permutation.all_permutations(n)
    ring = BoundedWCFactorAlgebra()

    for perm1 in perms:
        for perm2 in perms:
            elem1 = ring.full_groth_elem(perm1, n + 2, Gx._beta)
            elem2 = ring.full_groth_elem(perm2, n + 2, Gx._beta)
            prod = (elem1 * elem2).to_wc_graph_ring_element()
            expected_prod = Gx(perm1) * Gx(perm2)
            assert all(expected_prod.get(wc_result.perm, 0) == c * Gx._beta**(wc_result.perm.inv - perm1.inv - perm2.inv) for wc_result, c in prod.items() if wc_result.is_reduced), f"Error: Bounded WC factor algebra multiplication mismatch for permutations {perm1} and {perm2}, {prod=} {expected_prod=}"