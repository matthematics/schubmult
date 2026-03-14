from schubmult import *

if __name__ == "__main__":
    from schubmult.utils.schub_lib import all_grassmannian_rc_graphs
    import sys
    import itertools
    r = DualRCGraphRing()

    n = int(sys.argv[1])
    g_size = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    grasses = all_grassmannian_rc_graphs(n, g_size)
    
    perms = Permutation.all_permutations(n)

    #for perm in perms:
    for g1, g2, g3 in itertools.product(grasses, repeat=3):
        
        test1 = g1.left_squash(g2).left_squash(g3)
        test2 = g1.squash_product(g2).squash_product(g3)
        assert test1 == test2, f"Failure for g1={g1}, g2={g2}, g3={g3}\nLHS: {tuple(test1)}\nRHS: {tuple(test2)}"