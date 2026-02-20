from schubmult import *
from schubmult.utils.schub_lib import elem_sym_perms
from schubmult.combinatorial_reps.crystal_graph import CrystalGraphTensor
from sympy import pretty_print


P = PlacticAlgebra()

def unit_vector(i, n):
    vec = [0] * n
    vec[i - 1] = 1
    return tuple(vec)

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    #for perm1, perm2 in itertools.product(perms, repeat=2):
    for perm in perms:
        
        for k in range(1, n):
            for i in range(1, k + 1):
            # spickle = tuple(sorted([p for p, j in elem_sym_perms(perm, i, i) if j == i]))
            # spoikle = []
            # lw_fat = set()
#            for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n), RCGraph.all_rc_graphs(perm2, n)):
                for rc in RCGraph.all_rc_graphs(perm, n):
                    # elem_rc = RCGraph.principal_rc(~uncode([i]), n)
                    # elem_tab = elem_rc.hw_tab_rep()[1]
                    # tensor = CrystalGraphTensor(rc, next(iter(RCGraph.all_rc_graphs(uncode(unit_vector(k, k)), n, weight=unit_vector(i, n)))))
                    # # tab1 = rc1.hw_tab_rep()[1]
                    # # tab2 = rc2.hw_tab_rep()[1]
                    # hw, raise_seq = tensor.to_highest_weight()
                    # high_shape = [a for a in hw.crystal_weight if a != 0]
                    # yam_tab = Plactic.yamanouchi(high_shape).reverse_raise_seq(raise_seq)
                    # # assert tab1.rs_insert(*tab2.row_word) == yam_tab
                    # # print("yay")
                    # tab = rc.hw_tab_rep()[1]
                    # # tab2 = Plactic().rs_insert(i, *tab.row_word)
                    # cutcode = len(perm.trimcode)
                    # while cutcode > k:
                    #     low_part = rc.rowrange(0, cutcode - 1).normalize()
                    #     cutcode = len(low_part.trimcode)
                    # # low_part = low_part.resize(k)  
                    # monk_graph = next(iter(RCGraph.all_rc_graphs(uncode(unit_vector(k, k)), cutcode, weight=unit_vector(i, cutcode))))
                    # # else:
                    # #     monk_graph = next(iter(RCGraph.all_rc_graphs(uncode(unit_vector(len(low_part), len(low_part))), len(low_part), weight=unit_vector(i, len(low_part)))))
                    # pickle = low_part.squash_product(monk_graph)
                    # new_rc = RCGraph([*pickle[:k], *rc.rowrange(k, n).shiftup(k)])
                    fatbank = Sx(perm) * Sx(uncode(unit_vector(k, k)))
                    # while not new_rc.is_valid:# or new_rc.perm not in fatbank:

                    #     low_part = low_part.zero_out_last_row()
                    #     if len(low_part) < k:
                    #         low_part = low_part.resize(k)  
                    #         monk_graph = next(iter(RCGraph.all_rc_graphs(uncode(unit_vector(k, k)), k, weight=unit_vector(i, k))))
                    #     else:
                    #         monk_graph = next(iter(RCGraph.all_rc_graphs(uncode(unit_vector(len(low_part), len(low_part))), len(low_part), weight=unit_vector(i, len(low_part)))))
                    #     pickle = low_part.squash_product(monk_graph)
                    #     new_rc = RCGraph([*pickle[:k], *rc.rowrange(k, n).shiftup(k)])
                    
                    # tab2 = new_rc.hw_tab_rep()[1]
                    # #inserter = Plactic().rs_insert(i, *tab.row_word)
                    # #tab3_rc = rc.pieri_insert(k, [i])
                    # #assert tab2 == inserter
                    # pretty_print(new_rc)
                    if k < len(perm.trimcode):
                        #assert new_rc.is_valid, f"Failed on {perm} with k={k} and i={i}\nOriginal RC:\n{rc}\nOriginal tab:\n{tab}\nNew tab:\n{tab2}\nNew RC:\n{new_rc}\nNew RC tab:\n{new_rc.hw_tab_rep()[1]}"
                        low_part = rc.rowrange(0, k).normalize()
                        while len(low_part.perm.trimcode) > len(perm.trimcode):
                            low_part = low_part.zero_out_last_row()
                        cutlen = len(low_part)
                        cutcode = max(k, len(low_part.perm.trimcode))
                        low_part = low_part.resize(cutcode)
                        monk_graph = next(iter(RCGraph.all_rc_graphs(uncode(unit_vector(k, k)), cutcode, weight=unit_vector(i, cutcode))))
                        pickle = low_part.squash_product(monk_graph)
                        new_rc = RCGraph([*pickle[:min(cutlen,cutcode)], *rc.rowrange(min(cutlen,cutcode), n).shiftup(min(cutlen,cutcode))])
                        assert new_rc.is_valid, f"Failed on {perm} with"
                        assert new_rc.perm in fatbank
                    
                    #assert tab2 == tab3_rc.hw_tab_rep()[1], f"Failed on {perm} with k={k} and i={i}\nOriginal RC:\n{rc}\nOriginal tab:\n{tab}\nExpected tab:\n{tab2}\nNew RC:\n{tab3_rc}\nNew RC tab:\n{tab3_rc.hw_tab_rep()[1]}"
                # fitness = tensor.to_highest_weight()[0]
                # if fitness in lw_fat:
                #     continue
                # lw_fat.add(fitness)
                # the_bottom = rc.rowrange(0, i).normalize()
                # if len(the_bottom) < i:                    
                #     the_bottom = the_bottom.resize(i)
                # elem_rc = elem_rc.resize(i)
                # if len(the_bottom) > len(elem_rc):
                #     elem_rc = elem_rc.shiftup(len(the_bottom) - i).resize(len(the_bottom))
                # fartsingtab = next(iter(r(the_bottom) % r(elem_rc)))
                # fartsingtab = RCGraph([*fartsingtab[:i],*rc.rowrange(i, n).shiftup(i)]).hw_tab_rep()[1]
                #{rcc.hw_tab_rep()[1] for rcc in ((the_bottom)) % r(elem_extend)*r(rc.vertical_cut(i)[1])).keys()}
            #     inserter = tab1.rs_insert(*elem_tab.row_word)
            #     hw, raise_seq = tensor.to_highest_weight()
            #     tab2 = Plactic.yamanouchi([a for a in hw.crystal_weight if a != 0]).reverse_raise_seq(raise_seq)
            #     assert inserter == tab2, f"Failed on {perm} with i={i}"
            #     # assert inserter == fartsingtab, f"Failed on {perm} with i={i}\nRC:\n{rc}\nExpected tab:\n{inserter}\nActual tab:\n{fartsingtab}"
            #     # assert fartsingtab.hw_tab_rep()[1] == inserter, f"Failed on {perm} with i={i}\nRC:\n{rc}\nExpected tab:\n{inserter}\nActual tab:\n{fartsingtab.hw_tab_rep()[1]}"
            #     if uncode(tensor.crystal_weight) in spickle:
            #         spoikle.append(uncode(tensor.crystal_weight))
            #         assert inserter == RCGraph.principal_rc(uncode(tensor.crystal_weight)).hw_tab_rep()[1], f"Failed on {perm} with i={i}\nRC:\n{rc}\nExpected tab:\n{RCGraph.principal_rc(uncode(tensor.crystal_weight)).hw_tab_rep()[1]}\nActual tab:\n{elem_tab.rs_insert(*tab1.row_word)}"
            #         #assert RCGraph.principal_rc(uncode(fitness.crystal_weight)).hw_tab_rep()[1] == Plactic().rs_insert(*fitness.factors[1].hw_tab_rep()[1].row_word, *list(range(i, 0, -1))), f"Failed on {perm} with i={i}\nRC:\n{rc}\nFitness:\n{fitness}\nExpected tab:\n{RCGraph.principal_rc(uncode(fitness.crystal_weight)).hw_tab_rep()[1]}\nActual tab:\n{Plactic().rs_insert(*list(range(i, 0, -1)),*fitness.factors[1].hw_tab_rep()[1].row_word)}"
            # assert spickle == tuple(sorted(spoikle)), f"Failed on {perm} with i={i}\nExpected: {spickle}\nActual: {spoikle}"
            # new_rc = rc.monk_insert(i)
            # assert new_tab == new_rc.hw_tab_rep()[1], f"Failed on {perm}, row {i}\nOriginal RC:\n{rc}\nOriginal tab:\n{tab}\nNew tab:\n{new_tab}\nNew RC:\n{new_rc}\nNew RC tab:\n{new_rc.hw_tab_rep()[1]}"