from sympy import init_printing, pretty_print
from schubmult import *
from schubmult.rings.plactic_algebra import PlacticAlgebra
from schubmult.schub_lib.root_tableau import RootTableau

def rc_to_grass(rc):
    from schubmult.schub_lib.rc_graph import RCGraph
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
    T = ASx @ r
    g = GrassRCGraphRing()
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm):
            t_elem = ASx(rc.perm, n - 1)@ g(rc.grass)
            pretty_print(t_elem.coproduct())
    #Permutation.print_as_code = True
    # for perm1, perm2 in itertools.product(perms, repeat=2):
    #     prd1  = Sx(perm1)* Sx(perm2)
    #     hw_items = {}
    #     for w, coeff in prd1.items():
    #         for rc_w in RCGraph.all_hw_rcs(w, n - 1):
    #             for perm, coeff in prd1.items():
    #                 for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n - 1), RCGraph.all_rc_graphs(perm2, n - 1)):
    #                     squash_rc = rc1.squash_product(rc2)
    #                     if squash_rc.to_highest_weight()[0].length_vector == rc_w.length_vector and squash_rc.to_lowest_weight()[0].length_vector == rc_w.to_lowest_weight()[0].length_vector:
    #                         hw_items[w] = hw_items.get(w, r.zero) + coeff * r(squash_rc)
    #     for w in hw_items:
    #         print(f"{w}: {Sx(hw_items[w].polyvalue(Sx.genset))}")