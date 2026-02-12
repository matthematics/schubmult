from schubmult import *
from sympy import pretty_print

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        for rc in RCGraph.all_hw_rcs(perm, len(perm.trimcode)):
            chutes = rc.all_inverse_chute_moves()
            for chute in chutes:
                print("Inverse chute move in highest weight:")
                pretty_print(rc)
                print(chute)
                rc2 = rc.toggle_ref_at(*chute[0])
                rc2 = rc2.toggle_ref_at(*chute[1])
                pretty_print(rc2)
                assert rc2.to_highest_weight()[0] != rc.to_highest_weight()[0]
                print("IT IS A CRYSTAL BREAKER")
