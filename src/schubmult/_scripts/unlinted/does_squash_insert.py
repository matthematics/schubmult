from schubmult import *

P = PlacticAlgebra()

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    # cauchy
    # cauchy_prod = (Sx@r).one
    
    # for new_n in range(n - 1, 0, -1):
    #     new_cauchy_prod = (Sx@r).zero
    #     for deg in range(new_n + 1):
    #         if deg < new_n:
    #             elem_perm = uncode([0] * (new_n - deg) + [1] * deg)
    #             for rc in RCGraph.all_rc_graphs(elem_perm, new_n):
    #             cauchy_prod += sum(r(rc) for rc in RCGraph.all_rc_graphs(, new_n))
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, len(perm1.trimcode)), RCGraph.all_rc_graphs(perm2, len(perm2.trimcode))):
            # _, raise_seq1 = rc1.to_highest_weight()
            # _, raise_seq2 = rc2.to_highest_weight()           
            # raise_seq2 = [a + len(rc1) for a in raise_seq2]
            hw_rc1, tab1 = rc1.hw_tab_rep()
            hw_rc2, tab2 = rc2.hw_tab_rep()
            
            # the_prod = r(rc1)*r(rc2)
            # rat_prod = sum([coeff * (r@P)(rcc.hw_tab_rep()) for rcc, coeff in the_prod.items()])
            # all_hws = r(hw_rc1) * r(hw_rc2)
            # phanta_prod = (r@P).zero
            # for almost_hw_rc, coeff in all_hws.items():
            #     hw_rc, raise_seq = almost_hw_rc.to_highest_weight()
            #     high_shape = [a for a in hw_rc.length_vector if a != 0]
            #     yam_tab = Plactic.yamanouchi(high_shape).reverse_raise_seq(raise_seq).reverse_raise_seq(raise_seq1).reverse_raise_seq(raise_seq2)
            #     phanta_prod += coeff * (r@P)((hw_rc,yam_tab))
            # assert all(v == 0 for v in (rat_prod - phanta_prod).values()), f"Failed on {perm1} and {perm2}"
                

    #     grass_perms = [permg for permg in perms if len(permg.descents()) == 1 and max(permg.descents()) >= len(perm.trimcode) - 1]
    #     for rc in RCGraph.all_rc_graphs(perm, n - 1):
    #         tup1 = rc.hw_tab_rep()
    #         for grass_perm in grass_perms:
    #             for g_rc in RCGraph.all_rc_graphs(grass_perm, len(grass_perm.trimcode)):
    #                 tup2 = g_rc.hw_tab_rep()
    #                 squash_prod = next(iter((r(rc.resize(len(grass_perm.trimcode))) % r(g_rc)).keys()))
    #                 _, tab2 = squash_prod.hw_tab_rep()
    #                 assert tab2.row_word == tup2[1].rs_insert(*tup1[1].row_word).row_word, f"Failed on {perm} and {grass_perm}"
    # print("All squash   ")



    # does vex vex
        # for rc in RCGraph.all_hw_rcs(perm, len(perm.trimcode)):
        #     if not rc.perm.is_vexillary:
        #         assert len(rc[-1]) == 0
        #         while len(rc[-1]) == 0:
        #             rc = rc.zero_out_last_row()
        #         assert rc.perm.is_vexillary, f"Failed on {perm}\n{rc}"