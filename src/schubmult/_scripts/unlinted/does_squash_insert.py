from schubmult import *
from schubmult.schub_lib.crystal_graph import CrystalGraphTensor

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
        # if len(perm2.descents()) != 1 or len(perm2.trimcode) < len(perm1.trimcode):
        #     continue
        # if len(perm2.descents()) > 1 or len(perm2.trimcode) < len(perm1.trimcode):
        #     continue
        # if not perm1.is_dominant:
        #     continue
        #if any(perm1.has_pattern(pat) for pat in [[1,3,2]]) or len(perm1.descents()) > 1 or len(perm1.trimcode) >= len(perm2.trimcode):
        # if len(perm2.descents()) > 1:# or len(perm1.trimcode) >= len(perm2.trimcode):
        #     continue
        #if len(perm1.trimcode) > len(perm2.trimcode) or not perm2.is_vexillary:
        #if not perm2.is_vexillary or len(perm1.trimcode) > len(perm2.trimcode):
        #    continue
        print("Fat")
        # if len([a for a in perm1.trimcode if a != 0]) > 1:
        #     continue
        the_prod = Sx(perm1) * Sx(perm2)
        plactic_elem = P.zero
        for w, coeff in the_prod.items():
            plactic_elem += sum([coeff * P(rc.hw_tab_rep()[1]) for rc in RCGraph.all_rc_graphs(w, len(w.trimcode))])
        plactic_elem2 = P.zero
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n), RCGraph.all_rc_graphs(perm2, n)):
            # _, raise_seq1 = rc1.to_highest_weight()
            # _, raise_seq2 = rc2.to_highest_weight()           
            # raise_seq2 = [a + len(rc1) for a in raise_seq2]
            # hw_rc1, tab1 = rc1.hw_tab_rep()
            # hw_rc2, tab2 = rc2.hw_tab_rep()
            # plactic_elem2 += rc2.hw_tab_rep()[1].rs_insert(*rc1.hw_tab_rep()[1].row_word)
            tensor = CrystalGraphTensor(rc1, rc2)
            hw, raise_seq = tensor.to_highest_weight(length=min(len(perm1.trimcode), len(perm2.trimcode)))
            high_shape = [a for a in hw.crystal_weight if a != 0]
            yam_tab = Plactic.yamanouchi(high_shape).reverse_raise_seq(raise_seq)
            plactic_elem2 += P(yam_tab)
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
        assert all(v == 0 for v in (plactic_elem - plactic_elem2).values()), f"Failed on {perm1} and {perm2}\nplactic_elem: \n{plactic_elem}\nplactic_elem2: \n{plactic_elem2}"
    # for perm in perms:
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