from schubmult import *
from schubmult.rings.free_algebra import *

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.trimcode) for perm in perms]
    wr = WCGraphRing()
    GrothDual = FreeAlgebra(GrothendieckBasis)
    for perm1, perm2 in itertools.product(perms, repeat=2):
        for len1, len2 in itertools.product(range(perm1.max_descent, n), range(perm2.max_descent, n)):
                
            groth1 = GrothDual(perm1, len1)
            groth2 = GrothDual(perm2, len2)
            spinach1 = wr(WCGraph.principal_wc(perm1, len1))
            spinach2 = wr(WCGraph.principal_wc(perm2, len2))
            the_prod = spinach1 * spinach2
            real_prod = groth1 * groth2
            test_prod = 0
            for wc, v in the_prod.items():
                test_prod += v * GrothDual(wc.perm, len(wc))
            assert real_prod.almosteq(test_prod), f"Mismatch for {perm1=} {len1=} {perm2=} {len2=} {test_prod-real_prod=}\n{test_prod=}\n{real_prod=}"
            print(f"Success for {perm1=} {len1=} {perm2=} {len2=}")