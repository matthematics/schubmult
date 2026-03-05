from schubmult import *

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    r = RCGraphRing()
    perms = Permutation.all_permutations(n)
    for perm1, perm2 in itertools.permutations(perms, 2):
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n - 1), RCGraph.all_rc_graphs(perm2, n - 1)):
            cprd1 = r(rc1).to_free_algebra_element().coproduct()
            cprd2 = r(rc2).to_free_algebra_element().coproduct()
            cprdprd = cprd1 * cprd2
            try_cprd = (r(rc1)*r(rc2)).coproduct()
            try_cprd_free = sum([coeff * ASx(rc1a.perm, len(rc1a))@ASx(rc2a.perm, len(rc2a)) for (rc1a, rc2a), coeff in try_cprd.items()])
            assert try_cprd_free.almosteq(cprdprd), f"Failed for {rc1} and {rc2}, got {try_cprd_free-cprdprd}"
            print("Bacon potato")