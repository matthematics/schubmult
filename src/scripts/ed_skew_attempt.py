if __name__ == "__main__":
    import copy
    import sys
    from itertools import zip_longest

    import sympy
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
            crystals = {}
            highest_weights = set()
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
                
                for rc_w in RCGraph.all_rc_graphs(w, len(w.trimcode)):
                    pretty_print(rc_w)
                    # rc_hw = rc_w.to_highest_weight()[0]
                    # if rc_hw in crystals:
                    #     continue
                    if not rc_w.is_highest_weight:
                        continue
                    # if not rc_w.to_lowest_weight()[0].is_principal:
                    #     continue
                    high_weight = rc_w.length_vector
                    # if lw != tuple([a + b for a, b in zip_longest(prin_rc_u.length_vector, dom.rc_graph.length_vector, fillvalue=0)]):
                    #     print(f"{lw} != {tuple([a + b for a, b in zip_longest(prin_rc_u.length_vector, dom.rc_graph.length_vector, fillvalue=0)])}")
                    #     continue
                    w_tab = RootTableau.from_rc_graph(rc_w) # highest weight is the shape
                    
                    
                    skew_tab_set = set(NilPlactic.all_skew_ed_tableaux(w_tab.shape, u.antiperm, dom.shape))
                    for tb in skew_tab_set:
                        print('Found skew tab')
                        pretty_print(tb)
                        # the weight of the skew comes from w
                        assert tb.perm.inv == u.inv
                        # reduced_word = [len(u) - tb[box] for box in w_tab.iter_boxes_row_word_order if box in set(tb.iter_boxes)]
                        # new_grid = copy.deepcopy(w_tab._root_grid)
                        # index_list = []
                        # print(f"{w_tab=} {tb=} {u=}")
                        # for box in w_tab.iter_boxes:
                        #     if tb[box] == 0 :
                        #         new_grid[box] = None
                        #     else:
                        #         index_list.append(w_tab.order_grid[box])
                        # index_flatten = {b: len([a for a in index_list if a < b]) for b in index_list}
                        # for box in w_tab.iter_boxes:
                        #     if new_grid[box] is not None:
                        #         new_grid[box] = (u.right_root_at(index_flatten[w_tab.order_grid[box]], word=reduced_word), new_grid[box][1])
                        # compatible_seq = [w_tab[box][1] for box in tb.iter_boxes]
                        
                        # compatible_seq.sort()
                        # print(f"{reduced_word=} {compatible_seq=}")
                        # NEEEEEEEEEEEEEEEEEEEEEEEEEEEEEED DELETEL
                        u_tab = w_tab
                        box_grid = copy.deepcopy(u_tab._root_grid)
                        for box in tb.iter_boxes:
                            box_grid[box] = (u_tab[box][0], sympy.Integer(u_tab[box][1]))
                        u_tab_mover = RootTableau(box_grid)
                        u_roots = [u.right_root_at(index) for index in range(u.inv)]
                        while u_tab.perm.inv != u.inv:
                            for box in u_tab_mover.iter_boxes:
                                if not isinstance(u_tab_mover[box][1], sympy.Integer) and u_tab_mover[box][1] is not sympy.S.One:
                                    u_tab_test = u_tab.delete_box(box)
                                    if u_tab_test is not None:
                                        u_tab = u_tab_test
                                        u_tab_mover = u_tab_mover.delete_box(box)
                                        break
                        if u_tab.perm != u:
                            print("Skipping")
                            continue
                        pretty_print(u_tab)
                        pretty_print(tb.row_word)
                        print(f"Barfum {u.antiperm=}")
                        u_tab = u_tab.rectify()
                        u_hw_rc = u_tab.rc_graph
                        assert u_hw_rc.perm == u
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
                            tensor = CrystalGraphTensor(dom.rc_graph, u_tab2)
                            print(f"{tensor=}")
                            if not tensor.is_highest_weight:
                                print("Not highest weight")
                                continue
                            tc_elem = tensor
                            print(f"tc_elem=")
                            pretty_print(tc_elem)
                            if tc_elem in highest_weights:
                                print("Already there")
                                print(f"{highest_weights=}")
                                continue
                            pretty_print(dom.rc_graph)
                            assert tc_elem.crystal_weight == tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab2.length_vector, fillvalue=0)]), f"{tc_elem.crystal_weight=} vs {tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab2.length_vector, fillvalue=0)])}"
                            high_weight_check = tuple([a for a, b in zip_longest(high_weight, tc_elem.crystal_weight, fillvalue=0)])
                            low_weight_check = tuple([a for a, b in zip_longest(rc_w.to_lowest_weight()[0].length_vector, tc_elem.crystal_weight, fillvalue=0)])
                            if tc_elem.crystal_weight == high_weight_check and tc_elem.to_lowest_weight()[0].crystal_weight == low_weight_check:
                                crystals[rc_w] = crystals.get(rc_w, 0) + 1
                                print(f"{u=} {dom.perm=} {w=} {crystals=}")
                                highest_weights.add(tc_elem)
                            else:
                                print(f"{tc_elem.crystal_weight=}")
                                print(f"{high_weight_check=}")
                                # input()
                try:
                    any_cry = next(iter(crystals))
                    coeff = crystals[any_cry]
                    assert coeff == (Sx(dom.perm) * Sx(u)).get(w, 0), f"Fail at coeff check, {u=} {w=} {(Sx(dom.perm)*Sx(u)).get(w, 0)=} {coeff=}, {crystals=}"
                except AssertionError as e:
                    print(f"Failed coeff check at {u=} {w=} {dom.perm=}, expected {(Sx(dom.perm)*Sx(u)).get(w, 0)}, got {coeff}")
                    print("Crystals:")
                    for c in crystals:
                        pretty_print(c)
                    print("Highest weights:")
                    for c in highest_weights:
                        pretty_print(c)
                    raise
                print("Successful coeff check")
