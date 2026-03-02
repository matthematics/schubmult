from schubmult import *

if __name__ == "__main__":
    from schubmult.symbolic import expand, S
    from schubmult.abc import x
    from schubmult.utils.perm_utils import weak_compositions
    from schubmult.visualization import draw_pipe_dream_tikz
    from sympy import pretty_print, pretty, Mul
    import itertools
    import sys

    n = int(sys.argv[1])
    #m = int(sys.argv[2])
    r = RCGraphRing()
    #comps = weak_compositions(n, m)
    perms = Permutation.all_permutations(n)
    prod_dict = {}
    polystink = {}
    for perm in perms:
        # if perm1 != (1, 3, 4, 2) or perm2 != (4, 1, 3, 2):
        #     continue
        # if perm1.inv <= 1 or perm2.inv <= 1:
        #     continue
        #print(f"Computing product for {perm.trimcode}")
        inv_to_weight = {}
        perm_to_inv = {}
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            inv = rc.forest_invariant
            if rc.forest_weight not in inv_to_weight:
                inv_to_weight[rc.forest_weight] = inv
                perm_to_inv[inv] = rc.perm
            elif inv_to_weight[rc.forest_weight] != inv:
                continue
            polystink[rc.forest_weight] = polystink.get(rc.forest_weight, r.zero) + r(rc)
    for rcf1, rcf2 in itertools.product(polystink.keys(), repeat=2):
        print(f"Computing product for {rcf1} and {rcf2}")
        result = r.zero
        bigger = False
        for rc_bum1, rc_bum2 in itertools.product(polystink[rcf1].keys(), polystink[rcf2].keys()):
            perm1 = rc_bum1.perm
            perm2 = rc_bum2.perm
            cheat_prod = Sx(perm1) * Sx(perm2)
            stock_rc = {}
            stinkbat = {}
            for w in cheat_prod:
                if len(w) > n:
                    bigger = True
                for rc in RCGraph.all_rc_graphs(w, n - 1):
                    coprod = r(rc).coproduct()
                    
                    for (rc1, rc2), coeff in coprod.items():
                        # if rc1.perm not in stock_rc:
                        #     stock_rc[rc1.perm] = rc1
                        # # elif stock_rc[rc1.perm] != rc1:
                        # #     continue
                        # else:
                        #     continue
                        # if rc2.perm not in stock_rc:
                        #     stock_rc[rc2.perm] = rc2
                        # else:
                        #     continue

                        if rc1 == rc_bum1 and rc2 == rc_bum2:
                            #print("product of")
                            #print(rc1)
                            #print(draw_pipe_dream_tikz(rc1, max_size=len(rc1), flip_horizontal=True, top_labeled=False, show_refs=False, scale=1.0, outline_rows=None, clip_at_outline=True))
                            #print(rc2)
                            #print(draw_pipe_dream_tikz(rc2, max_size=len(rc2), flip_horizontal=True, top_labeled=False, show_refs=False, scale=1.0, outline_rows=None, clip_at_outline=True))
                            #print("is")
                            #print(coeff * r(rc))
                            #stinkbat[(rc1, rc2)] = stinkbat.get((rc1, rc2), r.zero) + coeff * r(rc)
                            result += coeff * r(rc)
                            #print(draw_pipe_dream_tikz(rc, max_size=len(rc), flip_horizontal=True, top_labeled=False, show_refs=False, scale=1.0, outline_rows=None, clip_at_outline=True))
        print(result)
        assert all(v >= 0 for v in result.values()), f"Negative coefficient found in product of {rcf1} and {rcf2}"
        if not bigger:
            assert expand(polystink[rcf1].polyvalue(Sx.genset) * polystink[rcf2].polyvalue(Sx.genset) - result.polyvalue(Sx.genset)) == S.Zero, f"Polynomial mismatch for product of {rcf1} and {rcf2}"
        print("yay paper")
        #for bacon, pig in stinkbat.items():
            # print("Fartulate")
            # print(bacon)
            # print("Grabonik")
            # print(pig)
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

