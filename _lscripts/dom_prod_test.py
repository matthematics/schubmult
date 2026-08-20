from schubmult import *

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = [perm for perm in Permutation.all_permutations(n) if (~perm).is_strict_dominant]
    for perm in perms:
        if perm.inv == 0:
            continue
        lambda1 = (~perm).trimcode
        sepdesc = [perm2 for perm2 in perms if perm2.inv > 0 and lambda1[-1] > (~perm2).trimcode[0]]
        for perm2 in sepdesc:
            the_schub = Sx(perm) * Sx(perm2)
            new_perm = next(iter(the_schub.keys()))
            new_part = [*lambda1, *(~perm2).trimcode]

            assert tuple(new_part) == tuple((~new_perm).trimcode), f"Failed for {perm} and {perm2}: got {new_perm} with code {(~new_perm).trimcode}, expected {new_part}"

            new_perm2 = perm * perm2
            assert new_perm == new_perm2, f"Failed for {perm} and {perm2}: got {new_perm} but expected {new_perm2}"
            print("Success")

