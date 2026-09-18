def all_reduced_subwords(reduced_word, u):
    if u.inv > len(reduced_word):
        return set()
    if u.inv == 0:
        return {()}
    ret_set = set()
    for index in range(len(reduced_word) - 1, -1, -1):
        a = reduced_word[index]
        if a - 1 in u.descents():
            new_u = u.swap(a - 1, a)
            old_set = all_reduced_subwords(reduced_word[:index], new_u)
            for subword in old_set:
                new_subword = (*subword, index)
                ret_set.add(new_subword)
    return ret_set

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

    class MarkedInteger(int):
        pass

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


                    # skew_tab_set = set(NilPlactic.all_skew_ed_tableaux(w_tab.shape, u.antiperm, dom.shape))

                        # w_tab = RootTableau.from_rc_graph(rc_w)
                    reduced_word = rc_w.reduced_word
                    for subword in all_reduced_subwords(reduced_word, u):
                        compatible_seq = [MarkedInteger(a) if index in subword else a for index, a in enumerate(rc_w.compatible_sequence)]
                        u_tab = RootTableau.root_insert_rsk(reduced_word, compatible_seq)
                        last_inv = 1000
                        while u_tab.perm.inv < last_inv:
                            last_inv = u_tab.perm.inv
                            for box in u_tab.iter_boxes:
                                if not isinstance(u_tab[box][1], MarkedInteger):
                                    u_tab_test = u_tab.delete_box(box)
                                    if u_tab_test is not None:
                                        u_tab = u_tab_test
                                        break
                        if u_tab.perm.inv > u.inv:
                            # didn't make it
                            continue
                        # if u_tab.perm != u:
                        #     print("Skipping")
                        #     print("U_tab = ")
                        #     pretty_print(u_tab)
                        #     print(f"{u=} {u_tab.perm=}")
                        #     continue
                        pretty_print(u_tab)
                        #pretty_print(tb.row_word)
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
                        hw_checked = set()
                        for u_tab2 in u_hw_rc.full_crystal:
                            tensor = CrystalGraphTensor(dom.rc_graph, u_tab2)
                            print(f"{tensor=}")
                            tc_elem = tensor.to_highest_weight()[0]
                            print("tc_elem=")
                            pretty_print(tc_elem)
                            if tc_elem in hw_checked:
                                print("Already checked")
                                print(f"{highest_weights=}")
                                continue
                            if tc_elem in highest_weights:
                                print("Already known highest weight mapped to some demazure crystal")
                                continue
                            u_tab_hw = tc_elem.factors[1]
                            hw_checked.add(tc_elem)
                            pretty_print(dom.rc_graph)
                            assert tc_elem.crystal_weight == tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab_hw.length_vector, fillvalue=0)]), f"{tc_elem.crystal_weight=} vs {tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab2.length_vector, fillvalue=0)])}"
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
