from schubmult import *
from schubmult.rings.free_algebra import *

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.trimcode) for perm in perms]
    wr = WCGraphRing()
    for comp1, comp2 in itertools.product(comps, repeat=2):
        grove1 = GroveDual(*comp1)
        grove2 = GroveDual(*comp2)
        spinach1 = wr.from_dict({k: v for k, v in wr.monomial(*comp1).items() if k.grove_weight == k.length_vector})
        spinach2 = wr.from_dict({k: v for k, v in wr.monomial(*comp2).items() if k.grove_weight == k.length_vector})
        the_prod = spinach1 * spinach2
        real_prod = grove1 * grove2
        test_prod = 0
        invariants = {}
        for wc, v in the_prod.items():
            if wc.grove_weight not in invariants:
                invariants[wc.grove_weight] = wc.grove_invariant
            elif invariants[wc.grove_weight] != wc.grove_invariant:
                continue
            test_prod += v * GroveDual(*wc.grove_weight)
        assert real_prod.almosteq(test_prod), f"Mismatch for {comp1=} {comp2=} {test_prod-real_prod=}\n{test_prod=}\n{real_prod=}"