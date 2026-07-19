from schubmult import *
from schubmult.rings.free_algebra import *
from schubmult.utils.tuple_utils import pad_tuple

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.trimcode) for perm in perms]
    wr = WCGraphRing()
    for comp01, comp02 in itertools.product(comps, repeat=2):
        for len1, len2 in itertools.product(range(len(comp01), n), range(len(comp02), n)):
            comp1 = pad_tuple(comp01, len1)
            comp2 = pad_tuple(comp02, len2)
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
            print(f"Success for {comp1=} {comp2=}")