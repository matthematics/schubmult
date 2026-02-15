from schubmult import *
from schubmult.utils.schub_lib import elem_sym_perms
from schubmult.schub_lib.crystal_graph import CrystalGraphTensor

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    for perm in perms:
        pp = perm.antiperm # skew ED
        up_perms = [up_perm for up_perm, diff in elem_sym_perms(perm, 1, 1) if diff == 1]
        left_tensor = RCGraph.principal_rc(uncode([1]), len(perm.trimcode) + 1)
        result = Sx.zero
        for up_perm in up_perms:
            up_rc = RCGraph.principal_rc(up_perm, len(perm.trimcode) + 1).to_highest_weight()[0]
            down_rc = up_rc.to_lowest_weight()[0]
            assert down_rc.is_principal
            
            principal_shape = NilPlactic.from_word(up_rc.perm_word).shape
            #print(principal_shape)
            tabs = NilPlactic.all_skew_ed_tableaux(principal_shape, perm, (1,))

            
            for tab in tabs:
                assert tab.is_increasing
                tab2 = tab.transpose().rectify().transpose()
                assert tab2.is_increasing
                #print(tab2)
                # bigot = tab2.hw_rc(len(perm.trimcode) + 1)
                # tab3 = NilPlactic.from_word([len(tab2.perm) - a for a in tab2.row_word])
                # assert tab3.perm == perm
                bigot = tab2.hw_rc(len(perm.trimcode) + 1)
                high_weights = set()    
                for rc in bigot.full_crystal:
                    tensor = CrystalGraphTensor(left_tensor, rc)
                    #hw_tens = tensor.to_highest_weight()[0]
                    
                    
                    lw_tens = tensor.to_lowest_weight()[0]

                    if tensor.crystal_weight == up_rc.length_vector and lw_tens.crystal_weight == down_rc.length_vector:
                        
                        print(f"Found match for {perm} in tensor product with {up_perm}")
                        print("Tab2")
                        print(tab2)
                        print("RC")
                        print(rc)
                        
                        if tensor in high_weights:
                            print(f"Already found this highest weight tensor before, skipping duplicate: {tensor}")
                            continue
                        high_weights.add(tensor)
                        
                        result += Sx(up_perm)
        assert result == Sx(perm) * Sx(uncode([1])), f"Result mismatch for {perm}: got {result}, expected {Sx(perm) * Sx(uncode([1]))}"
                    


