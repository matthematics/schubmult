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
            # needed?
            # if not u.bruhat_leq(dom.perm):
            #     continue
            def decompose_tensor_product(dom, u):
                crystals = {}
                highest_weights = set()
                perm_set = set((Sx(u)*Sx(dom.perm)).keys())
                for w in perm_set:
                    if not u.bruhat_leq(w):
                        continue
                    if not dom.perm.bruhat_leq(w):
                        continue

                    # print(f"Moving on to {u=} {w=} {dom.perm=}")
                    for rc_w in RCGraph.all_rc_graphs(w, len(w.trimcode)):
                        pretty_print(rc_w)
                        if not rc_w.is_highest_weight:
                            continue
                        high_weight = rc_w.length_vector
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

                            u_tab = u_tab.rectify()
                            u_hw_rc = u_tab.rc_graph
                            assert u_hw_rc.perm == u

                            hw_checked = set()
                            for u_tab2 in u_hw_rc.full_crystal:
                                tensor = CrystalGraphTensor(dom.rc_graph, u_tab2)
                                # print(f"{tensor=}")
                                tc_elem = tensor.to_highest_weight()[0]
                                pretty_print(tc_elem)
                                if tc_elem in hw_checked:
                                    # print("Already checked")
                                    # print(f"{highest_weights=}")
                                    continue
                                # needed!!!
                                if tc_elem in highest_weights:
                                    # print("Already known highest weight mapped to some demazure crystal")
                                    continue
                                u_tab_hw = tc_elem.factors[1]
                                # hw_checked.add(tc_elem)
                                pretty_print(dom.rc_graph)
                                assert tc_elem.crystal_weight == tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab_hw.length_vector, fillvalue=0)]), f"{tc_elem.crystal_weight=} vs {tuple([a + b for a,b in zip_longest(dom.rc_graph.length_vector, u_tab2.length_vector, fillvalue=0)])}"
                                high_weight_check = tuple([a for a, b in zip_longest(high_weight, tc_elem.crystal_weight, fillvalue=0)])
                                low_weight_check = tuple([a for a, b in zip_longest(rc_w.to_lowest_weight()[0].length_vector, tc_elem.crystal_weight, fillvalue=0)])
                                if tc_elem.crystal_weight == high_weight_check and tc_elem.to_lowest_weight()[0].crystal_weight == low_weight_check:
                                    crystals[rc_w] = crystals.get(rc_w, 0) + 1
                                    # print(f"{u=} {dom.perm=} {w=} {crystals=}")
                                    highest_weights.add(tc_elem)
                return crystals
            crystals = decompose_tensor_product(dom, u)

            try:
                test_prod = Sx(u)*Sx(dom.perm)
                for w in test_prod:
                    any_cry = None
                    for rcc in crystals:
                        if rcc.perm == w:
                            any_cry = rcc
                            break
                    coeff = crystals[any_cry]
                    assert coeff == test_prod[w], f"Fail at coeff check, {u=} {w=} {test_prod[w]=} {coeff=}, {crystals=}"
            except AssertionError as e:
                print(f"Failed coeff check at {u=} {w=} {dom.perm=}, expected {(Sx(dom.perm)*Sx(u)).get(w, 0)}, got {coeff}")
                print("Crystals:")
                for c in crystals:
                    pretty_print(c)
                # print("Highest weights:")
                # for c in highest_weights:
                #     pretty_print(c)
                raise
            print(f"Successful tensor decompose {dom.perm,u}")
