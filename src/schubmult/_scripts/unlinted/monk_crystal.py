
if __name__ == "__main__":
    import sys

    from schubmult import CrystalGraph, CrystalGraphTensor
    from schubmult import RCGraph
    import itertools
    from schubmult import *

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    for perm1, perm2, perm3 in itertools.product(perms, repeat=3):
        for rc1, rc2, rc3 in itertools.product(
            RCGraph.all_rc_graphs(perm1, n - 1),
            RCGraph.all_rc_graphs(perm2, n - 1),
            RCGraph.all_rc_graphs(perm3, n - 1),
        ):
            p1 = rc1.squash_product(rc2.squash_product(rc3))
            p2 = (rc1.squash_product(rc2)).squash_product(rc3)
            if p1 != p2:
                print("FAIL")
                print(perm1, perm2, perm3)
                print(rc1)
                print(rc2)
                print(rc3)
                print("Left:")
                print(p1)
                print("Right:")
                print(p2)
                sys.exit(1)
            print("PASS", perm1, perm2, perm3)