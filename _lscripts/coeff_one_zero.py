from schubmult import *
import itertools

def _t_it(t, perm):
    return uncode([a*t for a in perm.trimcode])

def main(n):
    perms = Permutation.all_permutations(n)
    for perm1, perm2 in itertools.product(perms, repeat=2):
        prd = Sx(perm1) * Sx(perm2)
        perm3 = next(iter([k for k, v in prd.items() if v == 1]))
        assert (Sx(_t_it(2, perm1))*Sx(_t_it(2, perm2))).get(_t_it(2,perm3),0) == 1

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    main(n)