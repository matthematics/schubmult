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