from schubmult import *

if __name__ == "__main__":
    from schubmult.symbolic import expand, S
    from schubmult.rings.combinatorial.dual_rc_graph_ring import DualRCGraphRing
    from schubmult.abc import x
    from schubmult.utils.perm_utils import weak_compositions
    from schubmult.visualization import draw_pipe_dream_tikz
    from sympy import pretty_print, pretty, Mul
    import itertools
    import sys

    n = int(sys.argv[1])
    #m = int(sys.argv[2])
    r = DualRCGraphRing()
    #comps = weak_compositions(n, m)
    perms = Permutation.all_permutations(n)
    prod_dict = {}
    for perm1, perm2 in itertools.product(perms, repeat=2):
        # if perm1 != (1, 3, 4, 2) or perm2 != (4, 1, 3, 2):
        #     continue
        if perm1.inv <= 1 or perm2.inv <= 1:
            continue
        print(f"Computing product for {perm1} and {perm2}")
        #cheat_prod = Sx(perm1) * Sx(perm2)
        result = r.zero
        stock_rc = {}
        stinkbat = {}
        #print(result)
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n - 1), RCGraph.all_rc_graphs(perm2, n - 1)):
            addbag = r(rc1) * r(rc2)
            result += addbag
            # pretty_print((rc1, rc2))
            # pretty_print(addbag)
        pretty_print(result)
        # for bacon, pig in stinkbat.items():
        #     pretty_print(bacon)
        #     pretty_print(pig)
            #print(draw_pipe_dream_tikz(bacon, max_size=len(bacon), flip_horizontal=True, top_labeled=False, show_refs=False, scale=1.0, outline_rows=None, clip_at_outline=True))
    # for comp in comps:
    #     if all(c == 0 for c in comp):
    #         continue
    #     pord = r.monomial(*comp)

    #     for rc, coeff in pord.items():
    #         cprd = r(rc).coproduct()
    #         for (rc1, rc2), coeff2 in cprd.items():                
    #             prod_dict[(rc1, rc2)] = prod_dict.get((rc1, rc2), r.zero) + coeff * coeff2 * r(rc)
    # poims = set()
    # for (rc1, rc2), coeff in prod_dict.items():
    #     if len(rc1.perm) > n + 1 or len(rc2.perm) > n + 1:
    #         continue
    #     if rc1.perm.inv <= 1 or rc2.perm.inv <= 1:
    #         continue
    #     if (rc1.perm, rc2.perm) not in poims:
            
    #         sm = r.zero
    #         for rc11, rc22 in itertools.product(RCGraph.all_rc_graphs(rc1.perm, n), RCGraph.all_rc_graphs(rc2.perm, n)):
    #             sm += prod_dict.get((rc11, rc22), r.zero)
    #         if all(len(rc.perm) <= n for rc in sm):
    #         #    print(f"Product for {rc1.perm} and {rc2.perm} is {sm}")
    #             print(f"Computing product for {rc1.perm} and {rc2.perm}")
    #             poims.add((rc1.perm, rc2.perm))
    #             pretty_print(sm)

