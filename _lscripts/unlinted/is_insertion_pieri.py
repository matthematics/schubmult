if __name__ == "__main__":
    from schubmult import *
    import sys
    from itertools import combinations_with_replacement

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    for k in range(n):
        
        for p in range(5):
            seqs = [tuple(reversed(seq)) for seq in combinations_with_replacement(range(1, k), p)]
            for perm in perms:
                all_rcs = {}
                for seq in seqs:
                    for rc in RCGraph.all_rc_graphs(perm, max(k, len(perm.trimcode))):
                        new_rc = rc.kogan_kumar_insert(k, seq)
                        if new_rc in all_rcs:
                            raise AssertionError(f"Error: Duplicate RC graph generated for permutation {perm} inserting at k={k} with seq={seq}:\n{new_rc}\nExisting seq: {all_rcs[new_rc][1]}\nExisting RC {all_rcs[new_rc][0]} \nNew seq: {seq}")
                        all_rcs[new_rc] = (rc, seq)