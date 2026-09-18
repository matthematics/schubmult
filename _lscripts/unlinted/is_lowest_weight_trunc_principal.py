from schubmult import *
from sympy import pretty_print

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        dct = {}
        for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)):
            # if dct.get(rc.to_highest_weight()[0], None) is None:
            #     dct[rc.to_highest_weight()[0]] = rc.to_lowest_weight()[0]
            # else:
            #     assert dct[rc.to_highest_weight()[0]] == rc.to_lowest_weight()[0], f"Failed on {perm}"
            # assert rc.to_highest_weight()[0].length_vector == tuple(sorted(rc.to_lowest_weight()[0].length_vector, reverse=True)), f"Failed on {perm}\n{rc}\nHighest weight: {rc.to_highest_weight()[0].length_vector}\nLowest weight: {rc.to_lowest_weight()[0].length_vector}"
            # inverse_chute_moves = rc.all_inverse_chute_moves()
            # if len(inverse_chute_moves) > 0:
            #     for inv_chute in inverse_chute_moves:
            #         r, c = inv_chute[1]
            #         max_r_chute_col = max(chute[0][1] for chute in inverse_chute_moves if chute[1][0] == r)
            #         the_chute = next(chute for chute in inverse_chute_moves if chute[1][0] == r and chute[0][1] == max_r_chute_col)
            #         new_rc = rc.toggle_ref_at(*the_chute[0])
            #         new_rc = new_rc.toggle_ref_at(*the_chute[1])
            #         assert new_rc.to_highest_weight()[0] != rc, f"Failed on {perm} with chute {the_chute}"
            #         print(new_rc.to_highest_weight()[0].length_vector)
            #         print(rc.length_vector)
            # working_rc = rc
            # if not working_rc.perm.is_vexillary:
            #     if len(working_rc[-1]) == 0:
            #         while len(working_rc[-1]) == 0:
            #             working_rc = working_rc.zero_out_last_row()
            #         assert working_rc.is_principal
            #     else:
            #         while not working_rc.perm.is_vexillary:
            #             working_rc = working_rc.shiftup(1).normalize().zero_out_last_row()
            #         assert working_rc.is_principal
            # if len(rc.all_inverse_chute_moves()) == 0:
            #     assert uncode(rc.to_lowest_weight()[0].length_vector).is_vexillary, f"Failed on {perm}\n{rc}\nHighest weight: {rc.to_highest_weight()[0].length_vector}\nLowest weight: {rc.to_lowest_weight()[0].length_vector}"
            # top_rc, raise_seq = rc.to_top_rc()
            # #assert top_rc.is_principal, f"Failed on {perm}\n{rc}\n{top_rc}"
            # assert len(top_rc.all_inverse_chute_moves()) == 0, f"Failed on {perm}\n{rc}\n{top_rc}"
            for i in range(1, len(rc)):
                pancakes = rc.chute_raise(i)
                if pancakes is not None:
                    if rc.raising_operator(i) is not None:
                    #     assert pancakes.to_highest_weight()[0] != rc.to_highest_weight()[0], f"Failed on {perm} with chute raise at row {i}\n{rc}\n{pancakes}"
                    # else:
                        assert pancakes == rc.raising_operator(i), f"Failed on {perm} with chute raise at row {i}\n{rc}\n{pancakes}\n{rc.raising_operator(i)}"
                # while not working_rc.perm.is_vexillary:
                #     if len(working_rc[-1]) == 0:
                #         working_rc = working_rc.zero_out_last_row()
                #     else:
                #         working_rc = working_rc.shiftup(1).normalize().zero_out_last_row()
                # assert rc.to_lowest_weight()[0].length_vector == working_rc.resize(len(rc)).to_lowest_weight()[0].length_vector, f"Failed on {perm}\n{working_rc}\nHighest weight: {rc.to_highest_weight()[0].length_vector}\nLowest weight: {working_rc.to_lowest_weight()[0].length_vector}"

            # if len(chute_moves) > 0:
            #     for chute in chute_moves:
            #         r, c = chute[0]
            #         max_r_chute_col = max(chute[0][1] for chute in chute_moves if chute[0][0] == r)
            #         the_chute = next(chute for chute in chute_moves if chute[0][0] == r and chute[0][1] == max_r_chute_col)
            #         new_rc = rc.toggle_ref_at(*the_chute[0])
            #         new_rc = new_rc.toggle_ref_at(*the_chute[1])
            #         assert new_rc.to_highest_weight()[0] != rc.to_highest_weight()[0], f"Failed on {perm} with chute {the_chute}\n{rc}\n{new_rc}"

                    
                    
    print("All lowest weight truncation principal")