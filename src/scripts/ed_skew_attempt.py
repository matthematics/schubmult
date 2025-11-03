if __name__ == "__main__":
    import copy
    import sys

    from sympy import pretty_print

    from schubmult import CrystalGraphTensor, NilPlactic, Permutation, RCGraph, RootTableau, Sx
    
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    dominant_tabs = set()
    for perm in perms:
        if perm != perm.minimal_dominant_above():
            continue
        dominant_tabs.update([RootTableau.from_rc_graph(rc) for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode))])

    for u in perms:
        if u.inv == 0:
            continue
        prin_rc_u = RootTableau.from_rc_graph(RCGraph.principal_rc(u, len(u.trimcode)))
        for w in perms:
            if not u.bruhat_leq(w):
                continue
            for dom in dominant_tabs:
                if not dom.perm.bruhat_leq(w):
                    continue
                if (Sx(u)*Sx(dom.perm)).get(w, 0) == 0:
                    continue
                coeff = 0
                for rc_w in RCGraph.all_rc_graphs(w, len(w.trimcode)):
                    w_tab = RootTableau.from_rc_graph(rc_w) # highest weight is the shape
                    low_weight = rc_w.to_lowest_weight()[0].length_vector

                    skew_tab_set = NilPlactic.all_skew_ed_tableaux(w_tab.shape, u.antiperm, dom.shape)
                    for tb in skew_tab_set:
                        # the weight of the skew comes from w
                        assert tb.perm.inv == u.inv
                        reduced_word = [len(u) - a for a in tb.row_word]
                        new_grid = copy.deepcopy(w_tab._root_grid)
                        index_list = []
                        print(f"{w_tab=} {tb=} {u=}")
                        for box in w_tab.iter_boxes():
                            if tb[box] == 0 :
                                new_grid[box] = None
                            else:
                                index_list.append(w_tab.order_grid[box])
                        index_flatten = {b: len([a for a in index_list if a < b]) for b in index_list}
                        for box in w_tab.iter_boxes():
                            if new_grid[box] is not None:
                                new_grid[box] = (u.right_root_at(index_flatten[w_tab.order_grid[box]], word=reduced_word), new_grid[box][1])
                        u_tab = RootTableau(new_grid).rectify()
                        pretty_print(u_tab)
                        assert u_tab.perm == u, "Rectified tableau has wrong permutation"
                        print("Recified successfully")
                        