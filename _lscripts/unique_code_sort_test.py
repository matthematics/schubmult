from schubmult import *
from schubmult.rings.polynomial_algebra import *

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        poly = Schub(perm, n - 1).change_basis(MonomialBasis)
        scode = tuple(sorted(perm.pad_code(n - 1), reverse=True))
        assert poly.get(scode, 0) == 1, f"{perm}"