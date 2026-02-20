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
    r = RCGraphRing()
    g = GrassRCGraphRing()
    T = r @ g
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1), RCGraph.all_rc_graphs(perm2)):
            # prod = T((rc1.to_highest_weight()[0],tableau_to_rc(RootTableau.from_rc_graph(rc1).weight_tableau, len(rc1)))) * T((rc2.to_highest_weight()[0],tableau_to_rc(RootTableau.from_rc_graph(rc2).weight_tableau, len(rc2))))
            # prod_update = T.zero
            # for (rc11, rc22), coeff in prod.items():
            #     if coeff == 0:
            #         continue
            #     # if len(rc22.perm.descents()) > 1:
            #     #     continue
            #     prod_update += coeff * T((rc11.to_highest_weight()[0], tableau_to_rc(RootTableau.from_rc_graph(rc22).weight_tableau, len(rc22))))
            # rcprod = r(rc1)*r(rc2)
            # tryprod = T.zero
            # for rc, coeff in rcprod.items():
            #     tryprod += coeff*T((rc.to_highest_weight()[0],tableau_to_rc(RootTableau.from_rc_graph(rc).weight_tableau, len(rc))))
            # assert all(v==0 for v in (prod - tryprod).values()), f"Failed for {perm1} and {perm2} with products {prod} and {tryprod} {prod - tryprod}"
            # prod = T.ext_multiply(((r(rc1)*r(rc2)).to_highest_weight()[0]),rc_to_grass(r(rc1)*r(rc2)))
            # # flatten1 = T.zero
            # # for (rc11, rc22), coeff in prod.items():
            # #     flatten1 += coeff * T((rc11.to_highest_weight()[0], rc22))
            # prod2 = T.ext_multiply(r(rc1.to_highest_weight()[0])*r(rc2.to_highest_weight()[0]),
            # rc_to_grass(r(rc1))*rc_to_grass(r(rc2)))
            # prod2 = T.from_dict({(rc11.to_highest_weight()[0], rc_to_grass(rc22)): coeff for (rc11, rc22),coeff in prod2.items() if coeff != 0})
            # assert all(v for v in (prod - prod2).values()), f"Failed for {perm1} and {perm2}"
            # prod = ((r(rc1.grass) * r(rc2.grass))).grass
            # prod2 = (r(rc1)*r(rc2)).grass
            # prod_comp1 = (r(rc1.to_highest_weight()[0])*r(rc2.to_highest_weight()[0])).to_highest_weight()[0]
            # prod_comp2 = (r(rc1.grass)*r(rc2.grass)).grass
            # prod = T.ext_multiply(prod_comp1, prod_comp2)

            # prod2_comp1 = (r(rc1)*r(rc2)).to_highest_weight()[0]
            # prod2_comp2 = (r(rc1)*r(rc2)).grass
            # prod2 = T.ext_multiply(prod2_comp1, prod2_comp2)
            # prod1 = T.from_dict({(rr1.to_highest_weight()[0], rr2): coeff for (rr1, rr2), coeff in (T.ext_multiply(r(rc1.to_highest_weight()[0]), g(rc2)).items()) if len(rr2.perm.descents()) <= 1})
            # prod2 = T.from_dict({(rr1.to_highest_weight()[0], rr2.grass): coeff for (rr1, rr2), coeff in (T.ext_multiply(r(rc1), r(rc2)).items())})
            
            # try:
            #     assert all(v == 0 for v in (prod1 - prod2).values()), f"Failed for {perm1} and {perm2} with products {prod1} and {prod2} {prod1 - prod2}"
            # except AssertionError as e:
            #     print(f"Failed for {perm1} and {perm2} with products ")
            #     pretty_print(prod1)
            #     pretty_print(prod2) 
            #     pretty_print(prod1 - prod2)
            #     raise e
            topple = g.coproduct_on_basis(rc1.grass) * g.coproduct_on_basis(rc2.grass)
            bitoons = g(rc1.grass) * g(rc2.grass)
            topple2 = (g@g).zero
            for rc, coeff in bitoons.items():
                topple2 += coeff * g.coproduct_on_basis(rc)
            assert all(v == 0 for v in (topple - topple2).values()), f"Failed coproduct compatibility for {perm1} and {perm2}"
    # for perm in perms:
    #     cd = perm.trimcode
    #     schubdonk = FA(*cd).change_basis(SchubertBasis)
    #     picko = r.from_free_algebra_element(schubdonk)
    #     pretty_print(picko)
    #     fat_donky = r.one
    #     for a in cd:
    #         fat_donky *= r(RCGraph.one_row(a))
    #     pretty_print(fat_donky)