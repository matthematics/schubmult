from sympy import init_printing, pretty_print
from schubmult import *
from schubmult.rings.combinatorial.plactic_algebra import PlacticAlgebra
from schubmult.combinatorics.root_tableau import RootTableau

def rc_to_grass(rc):
    from schubmult.combinatorics.rc_graph import RCGraph
    if isinstance(rc, RCGraphRingElement):
        ret = rc.ring.zero
        for rc0, coeff in rc.items():
            ret += coeff * rc.ring(rc_to_grass(RCGraph(rc0)))
        return ret
    hw, raise_seq = rc.to_highest_weight()
    code = [r for r in reversed(hw.length_vector)]
    perm = uncode(code)
    top_rc = next(iter(RCGraph.all_rc_graphs(perm, len(rc), weight = tuple(reversed(code)))))
    return top_rc.reverse_raise_seq(raise_seq)

if __name__ == "__main__":
    # test module functionality

    import itertools
    import sys

    
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    # sKEW DIV DIFF WEIGHT
    # is dual pieri Cauchy?
    r = BoundedRCFactorAlgebra()
    # g = GrassRCGraphRing()
    # T = r @ g
    # for perm1, perm2 in itertools.product(perms, repeat=2):
    #     for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1), RCGraph.all_rc_graphs(perm2)):
    #         # topple = g.coproduct_on_basis(rc1.grass) * g.coproduct_on_basis(rc2.grass)
            # bitoons = rc1.grass) * g(rc2.grass)
            # topple2 = (g@g).zero
            # for rc, coeff in bitoons.items():
            #     topple2 += coeff * g.coproduct_on_basis(rc)
            # assert all(v == 0 for v in (topple - topple2).values()), f"Failed coproduct compatibility for {perm1} and {perm2}"
    # for perm in perms:
    #     cd = perm.trimcode
    #     schubdonk = FA(*cd).change_basis(SchubertBasis)
    #     picko = r.from_free_algebra_element(schubdonk)
    #     pretty_print(picko)
    #     fat_donky = r.one
    #     for a in cd:
    #         fat_donky *= r(RCGraph.one_row(a))
    #     pretty_print(fat_donky)
    import random
    perms2 = [*perms]
    random.shuffle(perms2)
    perms3 = [*perms]
    random.shuffle(perms3)
    size = n
    d = DualRCGraphRing()
    bad_patterns = [[4,1,3,2],[1,4,3,2],[3,1,4,2]]
    avoids_bad = lambda p: all(not p.has_pattern(pat) for pat in bad_patterns)
    for perm1, perm2, perm3 in itertools.product(perms, perms2, perms3):
    #for perm1, perm2 in itertools.product(perms, perms2):
        # if 0 in (perm1.inv, perm2.inv, perm3.inv) or any(not avoids_bad(p) for p in (perm1, perm2, perm3)):
        #     continue
        for (rc1, cem_dict1), (rc2, cem_dict2), (rc3, cem_dict3) in itertools.product(RCGraph.full_CEM(perm1,n).items(), RCGraph.full_CEM(perm2,n).items(), RCGraph.full_CEM(perm3,n).items()):
        #for (rc1, cem_dict1), (rc2, cem_dict2) in itertools.product(RCGraph.full_CEM(perm1,n).items(), RCGraph.full_CEM(perm2,n).items()):
            elem1, elem2, elem3 = r.from_tensor_dict(cem_dict1, size), r.from_tensor_dict(cem_dict2, size), r.from_tensor_dict(cem_dict3, size)
            # elem1, elem2 = r.from_tensor_dict(cem_dict1, size), r.from_tensor_dict(cem_dict2, size)
            test1 = (elem1*elem2)*elem3
            test2 = elem1 * (elem2 * elem3)
            # test1 = elem1.dual_product(elem2).dual_product(elem3)#.to_rc_graph_ring_element()
            # test2 = elem1.dual_product(elem2.dual_product(elem3))#.to_rc_graph_ring_element()
            if not test1.almosteq(test2):
                print(f"FAIL: perm1={perm1}, perm2={perm2}, perm3={perm3}")
                print(f"  (elem1 * elem2) * elem3: {test1}")
                print(f"  elem1 * (elem2 * elem3): {test2}")
                sys.exit(1)
            print(f"PASS: perm1={perm1}, perm2={perm2}, perm3={perm3}")
            #pretty_print(test1)
        # #filt = lambda p: d.from_dict({k: v for k, v in p.items() if avoids_bad(k.perm)})
        # for rc1, rc2, rc3 in itertools.product(RCGraph.all_rc_graphs(perm1, n - 1), RCGraph.all_rc_graphs(perm2, n - 1), RCGraph.all_rc_graphs(perm3, n - 1)):
        #     test1 = d(rc1) * (d(rc2) * d(rc3))
        #     test2 = (d(rc1) * d(rc2)) * d(rc3)
        #     if not test1.almosteq(test2):
        #         print(f"FAIL: perm1={perm1}, perm2={perm2}, perm3={perm3}")
        #         print(f"  d(rc1) * (d(rc2) * d(rc3)): {test1}")
        #         print(f"  (d(rc1) * d(rc2)) * d(rc3): {test2}\n{rc1=}\n{rc2=}\n{rc3=}")
        #         sys.exit(1)