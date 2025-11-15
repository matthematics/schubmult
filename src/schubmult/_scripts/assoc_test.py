from sympy import init_printing, pretty_print



def is_isomorphic_crystal(tensor, rc_graph):
    hw1 = tensor.to_highest_weight()[0]
    hw2 = rc_graph.to_highest_weight()[0]
    if hw1.crystal_weight != hw2.crystal_weight:
        return False
    lw1 = hw1.to_lowest_weight()[0]
    lw2 = hw2.to_lowest_weight()[0]
    if lw1.crystal_weight != lw2.crystal_weight:
        return False
    return True


if __name__ == "__main__":
    # test module functionality

    import itertools
    import sys

    from schubmult.abc import x
    from symengine import S
    from sympy import pretty_print

    from schubmult import Permutation, RCGraph, RCGraphRing, RootTableau, Sx, CrystalGraphTensor
    # from schubmultutils.perm_utils import artin_sequences

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    # sKEW DIV DIFF WEIGHT
    # is dual pieri Cauchy?
    dom_perms = {perm.minimal_dominant_above() for perm in perms}
    for dom in dom_perms:
        for v in perms:
            prod = Sx(dom) * Sx(v)
            for w in prod:
                print(f"{w.trimcode}")
                coeff = 0
                for dom_rc in RCGraph.all_rc_graphs(dom, n-1):
                    for v_rc in RCGraph.all_rc_graphs(v, n-1):
                        dp_ret = v_rc.dualpieri(dom, w)
                        if len(dp_ret) > 0:
                            coeff += 1
                # for rc in RCGraph.all_hw_rcs(v, n-1):
                #     dom_rc = RCGraph.principal_rc(dom, n-1)
                #     tensor_hw = set()
                #     for v_rc in rc.full_crystal:
                #         tensor = CrystalGraphTensor(dom_rc, v_rc)
                #         tensor_hw.add(tensor.to_highest_weight()[0])
                    
                #     for thw in tensor_hw:
                #         for indiv_tensor in thw.full_crystal:
                #             dp_ret = indiv_tensor.factors[1].dualpieri(dom, w)
                #             if len(dp_ret) > 0:
                #                 # check if crystal isomorphic
                #                 w_rcs = RCGraph.all_rc_graphs(w, n-1, weight=indiv_tensor.crystal_weight)
                #                 assert len(w_rcs) > 0, f"No RC graphs for {w} with weight {indiv_tensor.crystal_weight}"
                #                 found_isomorphic = False
                #                 for w_rc in w_rcs:
                #                     if is_isomorphic_crystal(indiv_tensor, w_rc):
                #                         found_isomorphic = True
                #                         break
                #                 assert found_isomorphic, f"No isomorphic crystal found for tensor {indiv_tensor} and w {w}"
                #                 print("Successfully verified isomorphic crystal for tensor and w")
                #                 assert len(w_rcs) == 1
                #                 print("This crystal is unique")
                #                coeff += 1
                print(f"Coefficient for {dom} * {v} to {w} is {coeff}, should be {prod[w]}")
                assert coeff == prod[w], f"Coefficient mismatch for {dom} * {v} to {w}: got {coeff}, expected {prod[w]}"
                print("Verified")

            