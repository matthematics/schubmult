from schubmult import *

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n):
            print(rc.perm_word)
            lbs, dec = rc.omega_invariant
            assert lbs.is_valid, f"Error: LBS {lbs} from RC graph {rc} is not valid"
            assert dec.is_valid, f"Error: Dec {dec} from RC graph {rc} is not valid"
            word = [0] * perm.inv
            for node, label in dec.inorder_traversal:
                word[label - 1] = lbs(node.index).primary
            print(word)