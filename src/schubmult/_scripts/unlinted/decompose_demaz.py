from schubmult import *
from sympy import pretty_print

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    for perm1, perm2 in itertools.combinations(perms, 2):
        if perm1.is_dominant or perm2.is_dominant:
            continue
        for rc1, rc2 in itertools.product(RCGraph.all_hw_rcs(perm1, n - 1), RCGraph.all_hw_rcs(perm2, n - 1)):
            weight1, dperm1 = rc1.classify_demazure_crystal()
            weight2, dperm2 = rc2.classify_demazure_crystal()
            pretty_print(r(rc1) @ r(rc2))
            print(f"Demazure crystal of {perm1} with dominant weight {weight1} and lowest weight {dperm1}")
            print(f"Demazure crystal of {perm2} with dominant weight {weight2} and lowest weight {dperm2}")
            print()
            print(f"Decomposes? {Permutation.does_demazure_crystal_tensor_decompose(weight1, dperm1, weight2, dperm2)}")
            print(f"Other direction? {Permutation.does_demazure_crystal_tensor_decompose(weight2, dperm2, weight1, dperm1)}")