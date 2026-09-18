from schubmult import *
from schubmult.utils.schub_lib import elem_sym_perms
from schubmult.combinatorics.crystal_graph import CrystalGraphTensor

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    mulperm = uncode([1,1])
    tensperm = uncode([1,1])
    e = EGPlacticRing()
    for perm in perms:
        # pp = perm.antiperm # skew ED
        # up_perms = [up_perm for up_perm, diff in elem_sym_perms(perm, 1, len(mulperm.trimcode)) if diff == 1]
        #up_perms = [up_perm for up_perm, diff in elem_sym_perms(perm, mulperm.inv, len(mulperm.trimcode)) if diff == mulperm.inv]
        up_perms = tuple((Sx(perm) * Sx(mulperm)).keys())
        
        result = Sx.zero
        for up_perm in up_perms:
            graph_len = max(len(up_perm.trimcode), len(perm.trimcode) + 1) + 1
            left_tensor = RCGraph.principal_rc(tensperm, graph_len)
            down_rc = RCGraph.principal_rc(up_perm, graph_len)
            up_rc = down_rc.to_highest_weight()[0]
            assert down_rc.is_principal
            
            principal_shape = NilPlactic.from_word(up_rc.perm_word).shape
            principal_shape_t = NilPlactic.from_word(up_rc.perm_word).transpose().shape
            #print(principal_shape)
            tabs = NilPlactic.all_skew_ed_tableaux(principal_shape, perm, tuple((~mulperm).trimcode))
            used_set = set()
            if len(tabs) == 0:
                # print(f"No skew ED tableaux found for {perm} in tensor product with {up_perm}")
                print(f"No skew ED tableaux found for {perm} in tensor product with {up_perm}")
                continue
            for tab in tabs:
                assert tab.is_increasing
                # tab2 = tab.transpose().rectify().transpose()
                # assert tab2.is_increasing
                # #print(tab2)
                # # bigot = tab2.hw_rc(len(perm.trimcode) + 1)
                # # tab3 = NilPlactic.from_word([len(tab2.perm) - a for a in tab2.row_word])
                # # assert tab3.perm == perm
                # # bigot = tab2.hw_rc(len(perm.trimcode) + 1)
                #bigot = tab.hw_rc(graph_len)
                low_weight = down_rc.length_vector
                outer_shape, inner_shape = tab.transpose().skew_shape
                tab2 = tab.rectify()
                # if tab2 in used_set:
                #     print(f"Already rectified this tableau before, skipping duplicate: {tab2}")
                #     continue
                # used_set.add(tab2)
                the_ss_tabs = Plactic.all_ss_tableaux(principal_shape_t, len(perm.trimcode), mulperm.trimcode)
                for ss_tab in the_ss_tabs:
                    pokle = ss_tab.rectify()
                    # if pokle.shape != tab2.transpose().shape:
                    #     continue
                    # print("Tab:")
                    # print(tab)                    
                    # print("SS Tab:")
                    # print(ss_tab)                    
                    # print("Pokle:")
                    # print(pokle)
                    #if tuple([a + b for a, b in itertools.zip_longest(pokle.crystal_weight, left_tensor.length_vector, fillvalue=0)]) == low_weight:
                    try:
                        elem = e.key_to_rc_graph(((tab2, graph_len), pokle)).resize(graph_len)
                    except ValueError as E:
                        # print(E)
                        continue
                    except AssertionError as F:
                        # print(F)
                        continue
                    #if elem.is_lowest_weight:
                    
                    if elem.is_valid:
                        tens = CrystalGraphTensor(left_tensor, elem)
                        #assert tens.crystal_length() == graph_len, f"Tensor crystal length mismatch: expected {graph_len}, got {tens.crystal_length}, {tens=}"
                        tens_crystal_weight = tuple([a + b for a, b in itertools.zip_longest(pokle.crystal_weight, left_tensor.length_vector, fillvalue=0)])
                        assert len(tens_crystal_weight) == graph_len, f"Tensor crystal weight length mismatch: expected {graph_len}, got {len(tens_crystal_weight)}, {tens=}"
                        if pokle not in used_set and (all(left_tensor.phi(i) > pokle.epsilon(i) or elem.lowering_operator(i) is None for i in range(1, graph_len))) and tens_crystal_weight == low_weight:# and tens.to_highest_weight()[0].crystal_weight == up_rc.length_vector:
                        # if tens.crystal_weight == low_weight and tens.is_lowest_weight and tens not in used_set:# and tens.to_highest_weight()[0].crystal_weight == up_rc.length_vector:
                            result += Sx(up_perm)
                            used_set.add(pokle)
                # high_weights = set()    
                # for rc in bigot.full_crystal:
                #     tensor = CrystalGraphTensor(left_tensor, rc.resize(graph_len))
                #     #hw_tens = tensor.to_highest_weight()[0]
                    
                    
                #     hw_tens = tensor.to_highest_weight()[0]
                #     the_tenst_weight = [*tensor.crystal_weight]
                    
                #     the_hw_weight = [*hw_tens.crystal_weight]
                    
                #     assert len(tensor.crystal_weight) == graph_len, f"Tensor crystal weight length mismatch: expected {graph_len}, got {len(tensor.crystal_weight)}, {tensor=}"
                #     if tuple(the_tenst_weight) == down_rc.length_vector and tensor.is_lowest_weight and tuple(the_hw_weight) == up_rc.length_vector:
                        
                #         # print(f"Found match for {perm} in tensor product with {up_perm}")
                #         # # print("Tab2")
                #         # # print(tab)
                #         # print("RC")
                #         # print(rc)
                        
                #         if tensor in high_weights:
                #             print(f"Already found this highest weight tensor before, skipping duplicate: {tensor}")
                #             continue
                #         high_weights.add(tensor)
                        
                #         result += Sx(up_perm)
                    # the_tenst_weight[0] -= 1
                    # the_tenst_weight[1] += 1
                    # the_hw_weight[0] -= 1
                    # the_hw_weight[1] += 1
                    # if tuple(the_tenst_weight) == down_rc.length_vector and tuple(the_hw_weight) == up_rc.length_vector:
                        
                    #     print(f"Found match for {perm} in tensor product with {up_perm}")
                    #     # print("Tab2")
                    #     # print(tab)
                    #     print("RC")
                    #     print(rc)
                        
                    #     if tensor in high_weights:
                    #         print(f"Already found this highest weight tensor before, skipping duplicate: {tensor}")
                    #         continue
                    #     high_weights.add(tensor)
                        
                    #     result += Sx(up_perm)
        try:
            assert result == Sx(perm) * Sx(mulperm), f"Result mismatch for {perm}: got {result}, expected {Sx(perm) * Sx(mulperm)}"
            print("Success for", perm)
        except AssertionError as e:
            print(e)
            print("Perm:", perm)
            print("Up perms:", up_perms)
            for tab in tabs:
                print("Tab:")
                print(tab)
                print("Rect")
                print(tab.rectify())
                bigot = tab.hw_rc(graph_len).resize(graph_len)
                print("Bigot:")
                print(bigot)
            raise
                
                    


