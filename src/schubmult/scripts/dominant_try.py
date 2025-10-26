
from schubmult.rings.ck_ring import CoxeterKnuthRing


if __name__ == "__main__":
    from schubmult import *
    from schubmult.rings.rc_graph import RCGraph
    from schubmult.rings.crystal_graph import CrystalGraph, CrystalGraphTensor
    from schubmult.schub_lib.schub_lib import elem_sym_perms
    from sympy import pretty_print
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    dominant_graphs = {RCGraph.principal_rc(perm.minimal_dominant_above(), n-1).normalize() for perm in perms}

    
    for perm in perms:
        if perm.inv == 0:
            continue
        for dom in dominant_graphs:
            if dom.perm.inv == 0:
                continue
            if not dom.perm.bruhat_leq(perm):
                continue
            outer_shape = RCGraph.principal_rc(perm, len(perm.trimcode)).p_tableau.shape
            inner_shape = dom.weight_tableau.shape
            result = {}
            print(f"Trying {perm} {dom.perm}")
            tab_set = NilPlactic.all_skew_ed_tableaux(outer_shape, Permutation.w0(n), inner_shape)
            print("Got skew tableaux")
            tabs = {}
            for tab in tab_set:
                if tab.perm.inv != perm.inv - dom.perm.inv:
                    continue
                pretty_print("Skew tableau:")
                pretty_print(tab)
                rect_tab = tab.rectify()
                print("Rectified")
                pretty_print(rect_tab)
                ck_ring = CoxeterKnuthRing()
                if (Sx(dom.perm) * Sx(~rect_tab.perm)).get(perm, 0) and ((~rect_tab.perm) not in tabs or rect_tab in tabs[~rect_tab.perm]):
                    tabs[~rect_tab.perm] = tabs.get(~rect_tab.perm, set())
                    tabs[~rect_tab.perm].add(tab)
                    result[~rect_tab.perm] = result.get(~rect_tab.perm, 0) + 1
            print("Final result:")
            matches = {}
            print(f"{result=}")

            for k in perms:                
                product = (Sx(dom.perm) * Sx(k))
                if product.get(perm, 0) == result.get(k, 0):
                    matches[k] = True
                else:
                    matches[k] = False
                    print(f"Warning: mismatch! {k}: expected {product.get(perm, 0)}, got {result.get(k, 0)}")
                    print("Distinct tableaux:")
                    for tab in tabs.get(k, []):
                        pretty_print(tab)
                        pretty_print(tab.rectify())
                    input()
            print(f"Matches: {matches}")
            if any(not v for v in matches.values()):
                print("Mismatch found!")