if __name__ == "__main__":
    import copy
    import sys
    from itertools import zip_longest

    from sympy import pretty_print

    from schubmult import CrystalGraphTensor, NilPlactic, Permutation, RCGraph, RootTableau, Sx, uncode

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    dominant_tabs = set()
    for perm in perms:
        if perm != perm.minimal_dominant_above():
            continue
        dominant_tabs.update([RootTableau.from_rc_graph(rc) for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode))])

    for dom in dominant_tabs:
        for u in perms:
            if u.inv == 0:
                continue
            if not u.bruhat_leq(dom.perm):
                continue
            u_crystal_hw = RCGraph.principal_rc(u, len(u.trimcode)).to_highest_weight()[0]
            u_tab_crystal = RootTableau.from_rc_graph(u_crystal_hw)
            #prin_rc_u = RCGraph.principal_rc(u, len(u.trimcode)).to_highest_weight()[0]
            crystals = set()
            for w in perms:
                if not u.bruhat_leq(w):
                    continue

                if not dom.perm.bruhat_leq(w):
                    continue
                if (Sx(u)*Sx(dom.perm)).get(w, 0) == 0:
                    continue
                print(f"Moving on to {u=} {w=} {dom.perm=}")
                coeff = 0

                u_highest_weights = set()
                highest_weights = set()
                for rc_w in RCGraph.all_rc_graphs(w, len(w.trimcode)):
                    pretty_print(rc_w)
                    # rc_hw = rc_w.to_highest_weight()[0]
                    # if rc_hw in crystals:
                    #     continue
                    if not rc_w.is_highest_weight:
                        continue
                    if not rc_w.to_lowest_weight()[0].is_principal:
                        continue
                    high_weight = rc_w.length_vector
                    # if lw != tuple([a + b for a, b in zip_longest(prin_rc_u.length_vector, dom.rc_graph.length_vector, fillvalue=0)]):
                    #     print(f"{lw} != {tuple([a + b for a, b in zip_longest(prin_rc_u.length_vector, dom.rc_graph.length_vector, fillvalue=0)])}")
                    #     continue
                    w_tab = RootTableau.from_rc_graph(rc_w) # highest weight is the shape


                    skew_tab_set = set(NilPlactic.all_skew_ed_tableaux(w_tab.shape, u.antiperm, dom.shape))
                    for tb in skew_tab_set:
                        # the weight of the skew comes from w
                        assert tb.perm.inv == u.inv
                        reduced_word = [len(u) - a for a in tb.row_word]
                        # new_grid = copy.deepcopy(w_tab._root_grid)
                        # index_list = []
                        # print(f"{w_tab=} {tb=} {u=}")
                        # for box in w_tab.iter_boxes():
                        #     if tb[box] == 0 :
                        #         new_grid[box] = None
                        #     else:
                        #         index_list.append(w_tab.order_grid[box])
                        # index_flatten = {b: len([a for a in index_list if a < b]) for b in index_list}
                        # for box in w_tab.iter_boxes():
                        #     if new_grid[box] is not None:
                        #         new_grid[box] = (u.right_root_at(index_flatten[w_tab.order_grid[box]], word=reduced_word), new_grid[box][1])
                        # u_tab = RootTableau(new_grid)
                        pretty_print(tb)
                        pretty_print(tb.row_word)
                        print(f"Barfum {u.antiperm=}")
                        u_hw_rc = []
                        last_letter = -1
                        for r in reduced_word:
                            if r > last_letter:
                                u_hw_rc.append([])
                            u_hw_rc[-1].append(r)
                            last_letter = r
                        u_hw_rc = RCGraph([tuple(row) for row in u_hw_rc]).normalize()
                        print("Constructing tensor product")


                        # tensor_lw, _ = tensor.to_lowest_weight()
                        # lw_rc = tensor_lw.factors[1].rc_graph
                        # print("lw_rc=")
                        # pretty_print(tensor_lw)
                        # print(f"{tensor_lw.crystal_weight=}")
                        # print("Lowest weight of crystal")
                        # pretty_print(tensor_lw)
                        # print("low_weight")
                        # pretty_print(low_weight)

                        for u_tab2 in u_hw_rc.full_crystal:
                            tensor = CrystalGraphTensor(dom.rc_graph.resize(len(rc_w)), u_tab2.resize(len(rc_w)))
                            print(f"{tensor=}")
                            tc_elem = tensor.to_highest_weight()[0]
                            print(f"{tc_elem=}")
                            if (tc_elem, rc_w) in crystals:
                                continue
                            if tc_elem.crystal_weight == high_weight:
                                coeff += 1
                                crystals.add((tc_elem, rc_w))
                                print(f"{u=} {dom.perm=} {w=} {coeff=} {crystals=}")
                                highest_weights.add(tc_elem)
                            else:
                                print(f"{tc_elem.crystal_weight=}")
                                print(f"{high_weight=}")
                                # input()
                try:
                    assert coeff == (Sx(dom.perm) * Sx(u)).get(w, 0), f"Fail at coeff check, {u=} {w=} {(Sx(dom.perm)*Sx(u)).get(w, 0)=} {coeff=}, {crystals=}"
                except AssertionError as e:
                    print(f"Failed coeff check at {u=} {w=} {dom.perm=}, expected {(Sx(dom.perm)*Sx(u)).get(w, 0)}, got {coeff}")
                    print("Crystals:")
                    for c in crystals:
                        pretty_print(c)
                    raise
                print("Successful coeff check")
