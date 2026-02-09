from schubmult import *

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

    for perm in perms:
        weight_pair_set = set()
        for rc in RCGraph.all_hw_rcs(perm, len(perm.trimcode)):
            hw = rc.length_vector
            lw = rc.to_lowest_weight()[0].length_vector
            assert (hw, lw) not in weight_pair_set, f"Duplicate weight pair for {perm} with hw {hw} and lw {lw}"
            weight_pair_set.add((hw, lw))