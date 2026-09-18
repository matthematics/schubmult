if __name__ == "__main__":
    from schubmult import *
    import sys
    import numpy as np
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    for perm in perms:
        
        for monk_index in range(1, n):
            rc_length = len(perm.trimcode)
            if monk_index >= rc_length:
                continue
            for var_value in range(1, monk_index + 1):
                weight = [0] * rc_length
                weight[var_value - 1] = 1
                monk_rc = next(iter(RCGraph.all_rc_graphs(Permutation.ref_product(monk_index), rc_length, weight=weight)))
                for rc in RCGraph.all_rc_graphs(perm, rc_length):
                    first_rows = rc.vertical_cut(monk_index)[0].squash_product(monk_rc.rowrange(0, monk_index))
                    the_bpd = BPD.from_rc_graph(rc).resize(n)
                    first_rows_bpd = BPD.from_rc_graph(first_rows).resize(n)

                    

                    new_bpd = BPD(np.concatenate([first_rows_bpd.grid[:monk_index, :], the_bpd.grid[monk_index:, :]], axis=0))
                    print(new_bpd.to_rc_graph())
                    new_perm = new_bpd.to_rc_graph().perm
                    diff_indexes = [i for i in range(len(new_perm)) if new_perm[i] != perm[i]]
                    if len(diff_indexes) != 2 or diff_indexes[0] >= monk_index or diff_indexes[1] <= monk_index:
                        raise AssertionError(f"Error: Monk multiplication permutation mismatch for permutation {perm} with monk index {monk_index=} {new_perm=}")
