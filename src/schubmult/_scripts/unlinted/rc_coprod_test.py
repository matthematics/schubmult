from schubmult import *

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    r = RCGraphRing()
    perms = Permutation.all_permutations(n)
    for perm1, perm2 in itertools.permutations(perms, 2):
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n - 1), RCGraph.all_rc_graphs(perm2, n - 1)):
            cprd1 = r(rc1).coproduct()
            cprd2 = r(rc2).coproduct()
            cprdprd = cprd1 * cprd2
            try_cprd = (r(rc1)*r(rc2)).coproduct()
            assert try_cprd.almosteq(cprdprd), f"Failed for {rc1} and {rc2}, got {try_cprd-cprdprd}"
            print("Bacon potato")