
if __name__ == "__main__":
    from schubmult import *
    from schubmult.rings.rc_graph import RCGraph
    from schubmult.rings.crystal_graph import CrystalGraph, CrystalGraphTensor
    from schubmult.schub_lib.schub_lib import elem_sym_perms
    from sympy import pretty_print
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    dominant_graphs = {RCGraph.principal_rc(perm.minimal_dominant_above(), n-1) for perm in perms}

    
    for perm in perms:
        if perm.inv == 0:
            continue
        for dom in dominant_graphs:
            if dom.perm.inv == 0:
                continue
            inner_shape = dom.weight_tableau.shape
            perms2 = dict(Sx(perm)*Sx(dom.perm))
            perms2 = {k: v for k, v in perms2.items() if v != 0}
            print(f"Trying {perm} * {dom.perm}")
            for perm2 in perms2:
                print(f"Trying {perm2}")
                for rc2 in RCGraph.all_rc_graphs(perm2):
                    outer_shape = rc2.p_tableau.shape
                    print("Skew shape")
                    pretty_print(outer_shape)
                    print("--------------------")
                    pretty_print(inner_shape)
                    tab_set = NilPlactic.all_skew_ed_tableaux(outer_shape, perm, inner_shape)
                    print("Got skew tableaux")
                    for tab in tab_set:
                        pretty_print(tab)
                        rect_tab = tab.rectify()
                        print("Rectified")
                        pretty_print(rect_tab)
