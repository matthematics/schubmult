
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

    dominant_graphs = {RCGraph.principal_rc(perm.minimal_dominant_above(), n-1) for perm in perms}

    
    for perm in perms:
        if perm.inv == 0:
            continue
        for dom in dominant_graphs:
            if dom.perm.inv == 0:
                continue
            if not dom.perm.bruhat_leq(perm):
                continue
            if max(len(perm0) for perm0 in (Sx(dom.perm)*Sx(perm)).keys()) > n:
                continue
            result = {}
            for perm2 in perms:
                for rc2 in RCGraph.all_rc_graphs(perm2, n-1):
                    if rc2.is_dom_perm_yamanouchi(dom.perm, perm):
                        result[perm2] = result.get(perm2, 0) + 1
            
            matches = {}

            for k in perms:
                if k.inv == 0:
                    continue         
                product = (Sx(dom.perm) * Sx(k))
                if max(len(perm0) for perm0 in product.keys()) > n:
                    continue
                if product.get(perm, 0) == result.get(k, 0):
                    print("good")
                    matches[k] = True
                else:
                    matches[k] = False
                    print(f"Warning: mismatch! {k}: expected {product.get(perm, 0)}, got {result.get(k, 0)}")
                    print("Distinct rcs:")
                    for tab in tabs.get(k, []):
                        pretty_print(tab)
                        print(f"{(~tab.perm)=}")
                        print("Rectified:")
                        pretty_print(tab.rectify())
                        print(f"{~(tab.rectify().perm)=}")
                    print("Highest weight rcs:")
                    for tab in components.get(k, []):
                        print("hw_rc")
                        pretty_print(tab)
                        print("p_tableau")
                        pretty_print(tab.p_tableau)
                        print("q_tableau")
                        pretty_print(tab.q_tableau)
                        print("weight_tableau")
                        pretty_print(tab.weight_tableau)
                        print("highest_weight_tableau")
                        pretty_print(tab.weight_tableau.to_highest_weight(length=tab.crystal_length())[0])
                    for wt in w_tabs.get(k, []):
                        print("lw_rc")
                        pretty_print(wt)
                        print("p_tableau")
                        pretty_print(wt.p_tableau)
                        print("q_tableau")
                        pretty_print(wt.q_tableau)
                        print("weight_tableau")
                        pretty_print(wt.weight_tableau)
                        print("highest_weight_tableau")
                        pretty_print(wt.weight_tableau.to_highest_weight(length=wt.crystal_length())[0])
                        # pretty_print(tab)
                        # pretty_print(tab.rectify())
                    input()
                #print(f"Matches: {matches}")
                
                if any(not v for v in matches.values()):
                    print("Mismatch found!")